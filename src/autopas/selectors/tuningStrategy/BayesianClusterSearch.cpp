/**
 * @file BayesianClusterSearch.cpp
 * @author Fabio Gratl
 * @date 17.11.22
 */

#include "BayesianClusterSearch.h"

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"

autopas::BayesianClusterSearch::BayesianClusterSearch(const std::set<ContainerOption> &allowedContainerOptions,
                                                      const NumberSet<double> &allowedCellSizeFactors,
                                                      const std::set<TraversalOption> &allowedTraversalOptions,
                                                      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                                      const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                                      const std::set<Newton3Option> &allowedNewton3Options,
                                                      size_t maxEvidence, AcquisitionFunctionOption predAcqFunction,
                                                      const std::string &outputSuffix, size_t predNumLHSamples,
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
      _gaussianCluster({}, continuousDims, GaussianCluster::WeightFunction::evidenceMatchingScaledProbabilityGM, sigma,
                       _rng, GaussianCluster::defaultVecToString, outputSuffix),
      _neighbourFun([this](const Eigen::VectorXi &target) -> std::vector<std::pair<Eigen::VectorXi, double>> {
        return _encoder.clusterNeighboursManhattan1Container(target);
      }),
      _maxEvidence(maxEvidence),
      _predAcqFunction(predAcqFunction),
      _predNumLHSamples(predNumLHSamples),
      _firstTuningPhase(true),
      _currentIteration(0),
      _iterationScale(1. / (maxEvidence * iterationScalePerMaxEvidence)),
      _currentNumEvidence(0),
      _currentOptimalTime(std::numeric_limits<long>::max()),
      _fullSearch(allowedContainerOptions, {allowedCellSizeFactors.getMedian()}, allowedTraversalOptions,
                  allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options) {
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
    autopas::utils::ExceptionHandler::exception("BayesianClusterSearch: No valid configurations could be created.");
  }

  updateOptions();
  _gaussianCluster.setVectorToStringFun(
      [this](const GaussianModelTypes::VectorPairDiscreteContinuous &vec) -> std::string {
        return _encoder.convertFromCluster(vec).toString();
      });

  _currentConfig = _fullSearch.getCurrentConfiguration();
}

autopas::BayesianClusterSearch::~BayesianClusterSearch() = default;

const autopas::Configuration &autopas::BayesianClusterSearch::getCurrentConfiguration() const { return _currentConfig; }

void autopas::BayesianClusterSearch::addEvidence(long time, size_t iteration) {
  if (_firstTuningPhase) {
    _fullSearch.addEvidence(time, iteration);
  }

  // store optimal evidence
  if (time < _currentOptimalTime) {
    _currentOptimalTime = time;
    _currentOptimalConfig = _currentConfig;
  }

  // encoded vector
  auto vec = _encoder.convertToCluster(_currentConfig, iteration * _iterationScale);
  // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
  // represent a maximization problem
  _gaussianCluster.addEvidence(vec, -time * secondsPerMicroseconds);

  _currentIteration = iteration;
  ++_currentNumEvidence;
  _currentAcquisitions.clear();

  _traversalTimes[_currentConfig] = time;
}

long autopas::BayesianClusterSearch::getEvidence(Configuration configuration) const {
  return _traversalTimes.at(configuration);
}

void autopas::BayesianClusterSearch::reset(size_t iteration) {
  size_t iterationSinceLastEvidence = iteration - _currentIteration;
  if (iterationSinceLastEvidence * _iterationScale > suggestedMaxDistance) {
    AutoPasLog(warn, "BayesianClusterSearch: Time since the last evidence may be too long.");
  }

  _currentIteration = iteration;
  _currentNumEvidence = 0;
  _currentOptimalTime = std::numeric_limits<long>::max();
  _currentAcquisitions.clear();

  if (_firstTuningPhase) {
    _fullSearch.reset(iteration);
    _currentConfig = _fullSearch.getCurrentConfiguration();
  } else {
    tune();
  }
}

std::set<autopas::ContainerOption> autopas::BayesianClusterSearch::getAllowedContainerOptions() const {
  return _containerOptionsSet;
}

bool autopas::BayesianClusterSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _containerTraversalEstimatorOptions.size() == 1 and
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 1) and _dataLayoutOptions.size() == 1 and
         _newton3Options.size() == 1;
}

bool autopas::BayesianClusterSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerTraversalEstimatorOptions.empty() or
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or _dataLayoutOptions.empty() or
         _newton3Options.empty();
}

void autopas::BayesianClusterSearch::removeN3Option(Newton3Option badNewton3Option) {
  _fullSearch.removeN3Option(badNewton3Option);

  _newton3Options.erase(std::remove(_newton3Options.begin(), _newton3Options.end(), badNewton3Option),
                        _newton3Options.end());

  updateOptions();
  _currentAcquisitions.clear();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

void autopas::BayesianClusterSearch::updateOptions() {
  _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options,
                             *_cellSizeFactors);

  auto newRestrictions = _encoder.getDiscreteRestrictions();
  _gaussianCluster.setDimensions(std::vector<int>(newRestrictions.begin(), newRestrictions.end()));
}

bool autopas::BayesianClusterSearch::tune(bool currentInvalid) {
  if (currentInvalid) {
    _invalidConfigs.insert(_currentConfig);
  }

  // in the first tuning phase do a full search
  if (_firstTuningPhase) {
    if (_fullSearch.tune(currentInvalid)) {
      // continue with full-search
      _currentConfig = _fullSearch.getCurrentConfiguration();
      return true;
    } else {
      // continue with GaussianCluster
      _firstTuningPhase = false;
      _currentNumEvidence = 0;
    }
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = FeatureVector();
    return false;
  }

  // no more tunings steps
  if (_currentNumEvidence >= _maxEvidence) {
    // select best config of current tuning phase
    _currentConfig = _currentOptimalConfig;

    return false;
  }

  // try to sample a valid vector which is expected to yield a good acquisition
  for (size_t i = 0; i < maxAttempts; ++i) {
    if (_currentAcquisitions.empty()) {
      sampleAcquisitions(_predNumLHSamples, _predAcqFunction);
    }

    // test best vectors until empty
    while (not _currentAcquisitions.empty()) {
      auto best = _encoder.convertFromCluster(_currentAcquisitions.back());

      if (_invalidConfigs.find(best) == _invalidConfigs.end()) {
        // valid config found!
        _currentConfig = best;
        return true;
      } else {
        // invalid: dispose
        _currentAcquisitions.pop_back();
      }
    }

    // No valid configuration. This should rarely happen.
    AutoPasLog(debug, "Tuning could not generate a valid configuration.");
  }

  utils::ExceptionHandler::exception("BayesianClusterSearch: Failed to sample an valid FeatureVector");
  _currentConfig = FeatureVector();
  return false;
}

void autopas::BayesianClusterSearch::sampleAcquisitions(size_t n, AcquisitionFunctionOption af) {
  // create n lhs samples
  auto continuousSamples = _encoder.lhsSampleFeatureCluster(n, _rng, _currentIteration * _iterationScale);

  // calculate all acquisitions
  _currentAcquisitions = _gaussianCluster.sampleOrderedByAcquisition(af, _neighbourFun, continuousSamples);
}
