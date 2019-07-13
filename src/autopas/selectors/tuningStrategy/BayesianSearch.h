/**
 * @file BayesianSearch.h
 * @author Jan Nguyen
 * @date 12.06.19
 */

#pragma once

#include <limits>
#include <map>
#include <set>
#include <unordered_set>
#include "GaussianProcess.h"
#include "TuningStrategyInterface.h"
#include "autopas/selectors/ContainerSelector.h"
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
template <typename Particle, typename ParticleCell>
class BayesianSearch : public TuningStrategyInterface<Particle, ParticleCell> {
  /**
   * max limit of long as double
   */
  constexpr static double longMax = static_cast<double>(std::numeric_limits<long>::max());
  /**
   * The maximum number of attempts to sample an optimum.
   */
  constexpr static size_t maxAttempts = 10;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param predAcqFunction acquisition function used for prediction while tuning.
   * @param predNumSamples number of samples used for prediction while tuning.
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param lastAcqFunction acquisition function used for prediction of last tuning step.
   * @param lastNumSamples number of samples used for prediction of last tuning step.
   * @param seed seed of random number generator (should only be used for tests)
   */
  BayesianSearch(const std::set<ContainerOption> &allowedContainerOptions = allContainerOptions,
                 const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
                 const std::set<TraversalOption> &allowedTraversalOptions = allTraversalOptions,
                 const std::set<DataLayoutOption> &allowedDataLayoutOptions = allDataLayoutOptions,
                 const std::set<Newton3Option> &allowedNewton3Options = allNewton3Options, size_t maxEvidence = 10,
                 AcquisitionFunctionOption predAcqFunction = lcb, size_t predNumSamples = 1000,
                 AcquisitionFunctionOption lastAcqFunction = ucb, size_t lastNumSamples = 1000,
                 unsigned long seed = std::random_device()())
      : _containerOptions(allowedContainerOptions),
        _traversalOptions(allowedTraversalOptions),
        _dataLayoutOptions(allowedDataLayoutOptions),
        _newton3Options(allowedNewton3Options),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _traversalContainerMap(),
        _currentConfig(),
        _gp(1., std::vector<double>(FeatureVector::featureSpaceDims, 1.), 0.001),
        _invalidConfigs(),
        _maxEvidence(maxEvidence),
        _predAcqFunction(predAcqFunction),
        _predNumSamples(predNumSamples),
        _lastAcqFunction(lastAcqFunction),
        _lastNumSamples(lastNumSamples),
        _rng(seed) {
    /// @todo setting hyperparameters

    if (predNumSamples <= 0 or lastNumSamples <= 0) {
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

    tune();
  }

  inline const Configuration &getCurrentConfiguration() override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time) override {
    // transform time to a value between 0 and 1
    _gp.addEvidence(_currentConfig, time / longMax);
  }

  inline void reset() override {
    _gp.clear();
    tune();
  }

  inline bool tune(bool currentInvalid = false) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() override { return _containerOptions; }

  inline bool searchSpaceIsTrivial() override;

  inline bool searchSpaceIsEmpty() override;

  void addContainerSelector(/*const*/ ContainerSelector<Particle, ParticleCell> & /*containerSelector*/) override {
    AutoPasLog(trace, "ignoring containerSelector in BayesianSearch");
  }

 private:
  /**
   * Generate n samples and predict their corresponding
   * acquisition function. Return the FeatureVector which
   * generates the smallest value.
   * @param n numSamples
   * @param af acquisition function
   * @return optimal FeatureVector
   */
  inline FeatureVector sampleOptimalFeatureVector(size_t n, AcquisitionFunctionOption af);

  std::set<ContainerOption> _containerOptions;
  std::set<TraversalOption> _traversalOptions;
  std::set<DataLayoutOption> _dataLayoutOptions;
  std::set<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;

  std::map<TraversalOption, ContainerOption> _traversalContainerMap;

  FeatureVector _currentConfig;
  GaussianProcess<FeatureVector> _gp;
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;
  size_t _maxEvidence;
  AcquisitionFunctionOption _predAcqFunction;
  size_t _predNumSamples;
  AcquisitionFunctionOption _lastAcqFunction;
  size_t _lastNumSamples;
  Random _rng;
};

template <typename Particle, typename ParticleCell>
bool BayesianSearch<Particle, ParticleCell>::tune(bool currentInvalid) {
  if (currentInvalid) {
    _invalidConfigs.insert(_currentConfig);
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = FeatureVector();
    return false;
  }

  if (_gp.numEvidence() >= _maxEvidence) {
    // select predicted best config
    _currentConfig = sampleOptimalFeatureVector(_lastNumSamples, _lastAcqFunction);
    AutoPasLog(debug, "Selected Configuration {}", _currentConfig.toString());
    return false;
  }

  // select predicted best config for tuning
  _currentConfig = sampleOptimalFeatureVector(_predNumSamples, _predAcqFunction);
  return true;
}

template <typename Particle, typename ParticleCell>
FeatureVector BayesianSearch<Particle, ParticleCell>::sampleOptimalFeatureVector(size_t n, AcquisitionFunctionOption af) {
  for (size_t i = 0; i < maxAttempts; ++i) {
    // create n lhs samples
    std::vector<FeatureVector> samples = FeatureVector::lhsSampleFeatures(n, _rng, *_cellSizeFactors, _traversalOptions,
                                                                          _dataLayoutOptions, _newton3Options);

    // map container and calculate all acquisition function values
    std::map<FeatureVector, double> acquisitions;
    for (auto &sample : samples) {
      sample.container = _traversalContainerMap[sample.traversal];
      acquisitions[sample] = _gp.calcAcquisition(af, sample);
    }

    // sort by acquisition
    std::sort(samples.begin(), samples.end(), [&acquisitions](const FeatureVector &f1, const FeatureVector &f2) {
      return acquisitions[f1] < acquisitions[f2];
    });

    // find first valid configuration
    auto best = std::find_if(samples.begin(), samples.end(), [this](const FeatureVector &fv) {
      return _invalidConfigs.find(fv) == _invalidConfigs.end();
    });

    if (best != samples.end()) {
      // valid config found!
      return *best;
    } else {
      // No valid configuration. This should rarely happen.
      AutoPasLog(debug, "Tuning could not generate a valid configuration.");
    }
  }

  utils::ExceptionHandler::exception("BayesianSearch: Failed to sample an valid FeatureVector");
  return FeatureVector();
}

template <typename Particle, typename ParticleCell>
bool BayesianSearch<Particle, ParticleCell>::searchSpaceIsTrivial() {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _containerOptions.size() == 1 and (_cellSizeFactors->isFinite() && _cellSizeFactors->size() == 1) and
         _traversalOptions.size() == 1 and _dataLayoutOptions.size() == 1 and _newton3Options.size() == 1;
}

template <typename Particle, typename ParticleCell>
bool BayesianSearch<Particle, ParticleCell>::searchSpaceIsEmpty() {
  // if one enum is empty return true
  return _containerOptions.empty() or (_cellSizeFactors->isFinite() && _cellSizeFactors->size() == 0) or
         _traversalOptions.empty() or _dataLayoutOptions.empty() or _newton3Options.empty();
}

template <typename Particle, typename ParticleCell>
void BayesianSearch<Particle, ParticleCell>::removeN3Option(Newton3Option badNewton3Option) {
  _newton3Options.erase(badNewton3Option);

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

}  // namespace autopas
