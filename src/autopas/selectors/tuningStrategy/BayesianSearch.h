/**
 * @file BayesianSearch.h
 * @author Jan Nguyen
 * @date 12.06.2019
 */

#pragma once

#include <limits>
#include <map>
#include <set>
#include <unordered_set>

#include "GaussianModel/GaussianProcess.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/FeatureVector.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * Assume that the stochastic distribution of the execution time corresponds
 * to a Gaussian Process. This allows to estimate the 'gain' of testing a given
 * feature next.
 */
class BayesianSearch : public TuningStrategyInterface {
  /**
   * The maximum number of attempts to sample an optimum.
   */
  constexpr static size_t maxAttempts = 10;
  /**
   *
   */
  constexpr static double secondsPerMicroseconds = 1. / 1000000.;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param predAcqFunction acquisition function used for prediction while tuning.
   * @param predNumLHSamples number of samples used for prediction while tuning.
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param seed seed of random number generator (should only be used for tests)
   */
  explicit BayesianSearch(
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(), size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      size_t predNumLHSamples = 1000, unsigned long seed = std::random_device()())
      : _containerOptions(allowedContainerOptions),
        _traversalOptions(allowedTraversalOptions),
        _dataLayoutOptions(allowedDataLayoutOptions),
        _newton3Options(allowedNewton3Options),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _traversalContainerMap(),
        _currentConfig(),
        _invalidConfigs(),
        _rng(seed),
        _gaussianProcess(FeatureVector::oneHotDims, 0.01, _rng),
        _maxEvidence(maxEvidence),
        _predAcqFunction(predAcqFunction),
        _predNumLHSamples(predNumLHSamples) {
    if (predNumLHSamples <= 0) {
      utils::ExceptionHandler::exception(
          "BayesianSearch: Number of samples used for predictions must be greater than 0!");
    }

    // map each travesal to a container
    for (auto it = _traversalOptions.begin(); it != _traversalOptions.end();) {
      auto compatibleContainerOptions = compatibleTraversals::allCompatibleContainers(*it);
      std::set<ContainerOption> allowedAndCompatible;
      std::set_intersection(allowedContainerOptions.begin(), allowedContainerOptions.end(),
                            compatibleContainerOptions.begin(), compatibleContainerOptions.end(),
                            std::inserter(allowedAndCompatible, allowedAndCompatible.begin()));
      if (allowedAndCompatible.empty()) {
        // no compatible container found
        it = _traversalOptions.erase(it);
      } else {
        _traversalContainerMap[*it] = *allowedAndCompatible.begin();
        ++it;
      }
    }

    if (searchSpaceIsEmpty()) {
      autopas::utils::ExceptionHandler::exception("BayesianSearch: No valid configurations could be created.");
    }

    tune();
  }

  inline const Configuration &getCurrentConfiguration() const override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t iteration) override {
    // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
    // represent a maximization problem
    _gaussianProcess.addEvidence(_currentConfig.oneHotEncode(), -static_cast<double>(time) * secondsPerMicroseconds,
                                 true);
    _currentSamples.clear();
  }

  inline long getEvidence(Configuration configuration) const override {
    return -static_cast<long>(_gaussianProcess.predictMean(FeatureVector(configuration).oneHotEncode())
                              / secondsPerMicroseconds);
  }

  inline void reset(size_t iteration) override {
    _gaussianProcess.clear();
    _currentSamples.clear();
    tune();
  }

  inline bool tune(bool currentInvalid = false) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override { return _containerOptions; }

  inline bool searchSpaceIsTrivial() const override;

  inline bool searchSpaceIsEmpty() const override;

 private:
  /**
   * Generate n samples and sort them depending on
   * given acquistion function.
   * @param n numSamples
   * @param af acquisition function
   */
  inline void sampleAcquisitions(size_t n, AcquisitionFunctionOption af);

  std::set<ContainerOption> _containerOptions;
  std::set<TraversalOption> _traversalOptions;
  std::set<DataLayoutOption> _dataLayoutOptions;
  std::set<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;

  std::map<TraversalOption, ContainerOption> _traversalContainerMap;

  FeatureVector _currentConfig;
  std::vector<FeatureVector> _currentSamples;
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;

  Random _rng;
  GaussianProcess _gaussianProcess;
  size_t _maxEvidence;
  AcquisitionFunctionOption _predAcqFunction;
  size_t _predNumLHSamples;
};

bool BayesianSearch::tune(bool currentInvalid) {
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
    _currentConfig = FeatureVector::oneHotDecode(_gaussianProcess.getEvidenceMax());
    _currentConfig.container = _traversalContainerMap[_currentConfig.traversal];
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

void BayesianSearch::sampleAcquisitions(size_t n, AcquisitionFunctionOption af) {
  // create n lhs samples
  _currentSamples = FeatureVector::lhsSampleFeatures(n, _rng, *_cellSizeFactors, _traversalOptions, _dataLayoutOptions,
                                                     _newton3Options);

  // map container and calculate all acquisition function values
  std::map<FeatureVector, double> acquisitions;
  for (auto &sample : _currentSamples) {
    sample.container = _traversalContainerMap[sample.traversal];
    acquisitions[sample] = _gaussianProcess.calcAcquisition(af, sample.oneHotEncode());
  }

  // sort by acquisition
  std::sort(_currentSamples.begin(), _currentSamples.end(),
            [&acquisitions](const FeatureVector &f1, const FeatureVector &f2) {
              return acquisitions[f1] < acquisitions[f2];
            });
}

bool BayesianSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _containerOptions.size() == 1 and (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 1) and
         _traversalOptions.size() == 1 and _dataLayoutOptions.size() == 1 and _newton3Options.size() == 1;
}

bool BayesianSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerOptions.empty() or (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or
         _traversalOptions.empty() or _dataLayoutOptions.empty() or _newton3Options.empty();
}

void BayesianSearch::removeN3Option(Newton3Option badNewton3Option) {
  _newton3Options.erase(badNewton3Option);
  _currentSamples.clear();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

}  // namespace autopas
