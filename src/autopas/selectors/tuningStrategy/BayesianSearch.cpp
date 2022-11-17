/**
 * @file BayesianSearch.cpp
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "BayesianSearch.h"

#include "autopas/utils/StringUtils.h"

autopas::BayesianSearch::BayesianSearch(const std::set<ContainerOption> &allowedContainerOptions,
                                        const autopas::NumberSet<double> &allowedCellSizeFactors,
                                        const std::set<TraversalOption> &allowedTraversalOptions,
                                        const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                        const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                        const std::set<Newton3Option> &allowedNewton3Options, size_t maxEvidence,
                                        autopas::AcquisitionFunctionOption predAcqFunction, size_t predNumLHSamples,
                                        unsigned long seed)
    : _containerOptionsSet(allowedContainerOptions),
      _dataLayoutOptions(allowedDataLayoutOptions.begin(), allowedDataLayoutOptions.end()),
      _newton3Options(allowedNewton3Options.begin(), allowedNewton3Options.end()),
      _cellSizeFactors(allowedCellSizeFactors.clone()),
      _encoder(),
      _currentConfig(),
      _invalidConfigs(),
      _traversalTimes(),
      _rng(seed),
      _gaussianProcess(0, 0.01, _rng),
      _maxEvidence(maxEvidence),
      _predAcqFunction(predAcqFunction),
      _predNumLHSamples(predNumLHSamples) {
  if (predNumLHSamples <= 0) {
    utils::ExceptionHandler::exception(
        "BayesianSearch: Number of samples used for predictions must be greater than 0!");
  }

  for (const auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (const auto &traversalOption : allowedAndApplicable) {
      // if load estimators are not applicable LoadEstimatorOption::none is returned.
      const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
          loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
      for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
        _containerTraversalEstimatorOptions.emplace_back(containerOption, traversalOption, loadEstimatorOption);
      }
    }
  }

  if (searchSpaceIsEmpty()) {
    autopas::utils::ExceptionHandler::exception("BayesianSearch: No valid configurations could be created.");
  }

  _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options,
                             *_cellSizeFactors);
  _gaussianProcess.setDimension(_encoder.getOneHotDims());

  tune();
}

bool autopas::BayesianSearch::tune(bool currentInvalid) {
  if (currentInvalid) {
    _invalidConfigs.insert(_currentConfig);
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = FeatureVector();
    return false;
  }

  if (_gaussianProcess.numEvidence() >= _maxEvidence) {
    // select best config
    _currentConfig = _encoder.oneHotDecode(_gaussianProcess.getEvidenceMax());
    AutoPasLog(debug, "Selected Configuration {}", _currentConfig.toString());
    return false;
  }

  // try to sample a valid vector which is expected to yield a good acquisition
  for (size_t i = 0; i < maxAttempts; ++i) {
    if (_currentSamples.empty()) {
      sampleAcquisitions(_predNumLHSamples, _predAcqFunction);
    }

    // test best vectors until empty
    while (not _currentSamples.empty()) {
      if (_invalidConfigs.find(_currentSamples.back()) == _invalidConfigs.end()) {
        // valid config found!
        _currentConfig = _currentSamples.back();
        return true;
      } else {
        // invalid: dispose
        _currentSamples.pop_back();
      }
    }

    // No valid configuration. This should rarely happen.
    AutoPasLog(debug, "Tuning could not generate a valid configuration.");
  }

  utils::ExceptionHandler::exception("BayesianSearch: Failed to sample an valid FeatureVector");
  _currentConfig = FeatureVector();
  return false;
}

void autopas::BayesianSearch::sampleAcquisitions(size_t n, AcquisitionFunctionOption af) {
  // create n lhs samples
  _currentSamples = _encoder.lhsSampleFeatures(n, _rng);

  // map container and calculate all acquisition function values
  std::map<FeatureVector, double> acquisitions;
  for (auto &sample : _currentSamples) {
    acquisitions[sample] = _gaussianProcess.calcAcquisition(af, _encoder.oneHotEncode(sample));
  }

  // sort by acquisition
  std::sort(_currentSamples.begin(), _currentSamples.end(),
            [&acquisitions](const FeatureVector &f1, const FeatureVector &f2) {
              return acquisitions[f1] < acquisitions[f2];
            });
}

bool autopas::BayesianSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _containerTraversalEstimatorOptions.size() == 1 and
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 1) and _dataLayoutOptions.size() == 1 and
         _newton3Options.size() == 1;
}

bool autopas::BayesianSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerTraversalEstimatorOptions.empty() or
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or _dataLayoutOptions.empty() or
         _newton3Options.empty();
}

void autopas::BayesianSearch::removeN3Option(Newton3Option badNewton3Option) {
  _newton3Options.erase(std::remove(_newton3Options.begin(), _newton3Options.end(), badNewton3Option),
                        _newton3Options.end());
  _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options,
                             *_cellSizeFactors);

  _gaussianProcess.setDimension(_encoder.getOneHotDims());
  _currentSamples.clear();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

const autopas::Configuration &autopas::BayesianSearch::getCurrentConfiguration() const { return _currentConfig; }

void autopas::BayesianSearch::addEvidence(long time, size_t iteration) {
  // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
  // represent a maximization problem
  _gaussianProcess.addEvidence(_encoder.oneHotEncode(_currentConfig), -time * secondsPerMicroseconds, true);
  _currentSamples.clear();
  _traversalTimes[_currentConfig] = time;
}

long autopas::BayesianSearch::getEvidence(autopas::Configuration configuration) const {
  return _traversalTimes.at(configuration);
}

void autopas::BayesianSearch::reset(size_t iteration) {
  _gaussianProcess.clear();
  _currentSamples.clear();
  _traversalTimes.clear();
  tune();
}

std::set<autopas::ContainerOption> autopas::BayesianSearch::getAllowedContainerOptions() const {
  return _containerOptionsSet;
}

bool autopas::BayesianSearch::smoothedHomogeneityAndMaxDensityNeeded() const { return false; }
