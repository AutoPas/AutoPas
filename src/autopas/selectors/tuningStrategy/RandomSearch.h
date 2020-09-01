/**
 * @file RandomSearch.h
 * @author Jan Nguyen
 * @date 10.07.2019
 */

#pragma once

#include <map>
#include <set>

#include "TuningStrategyInterface.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Randomly test a given number of configurations and select the fastest.
 */
class RandomSearch : public TuningStrategyInterface {
 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param seed seed of random number generator (should only be used for tests)
   */
  RandomSearch(const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
               const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1. / 3., 2.),
               const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
               const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
               const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
               const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(),
               size_t maxEvidence = 10, unsigned long seed = std::random_device()())
      : _containerOptions(allowedContainerOptions),
        _traversalOptions(allowedTraversalOptions),
        _loadEstimatorOptions(allowedLoadEstimatorOptions),
        _dataLayoutOptions(allowedDataLayoutOptions),
        _newton3Options(allowedNewton3Options),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _currentConfig(),
        _invalidConfigs(),
        _rng(seed),
        _maxEvidence(maxEvidence),
        _searchSpaceSizeNoCSF(0),
        _numTestedConfigsNoCSF(0) {
    if (searchSpaceIsEmpty()) {
      autopas::utils::ExceptionHandler::exception("RandomSearch: No valid configurations could be created.");
    }

    // Allow invalid container-traversal-combinations, because so does tune()
    for (auto container : _containerOptions) {
      for (auto traversal : _traversalOptions) {
        auto applicableLoadEstimators =
            loadEstimators::getApplicableLoadEstimators(container, traversal, _loadEstimatorOptions);
        _searchSpaceSizeNoCSF += applicableLoadEstimators.size();
      }
    }
    _searchSpaceSizeNoCSF *= _dataLayoutOptions.size() * _newton3Options.size();

    // if fewer than _maxEvidence configurations exist in total, search through all of them and then stop.
    _maxEvidence = std::min(_searchSpaceSizeNoCSF, _maxEvidence);
    tune();
  }

  inline const Configuration &getCurrentConfiguration() const override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t iteration) override {
    Configuration testingConfig = Configuration(_currentConfig);
    testingConfig.cellSizeFactor = 1;
    if (_traversalTimes.find(testingConfig) == _traversalTimes.end()) {
      ++_numTestedConfigsNoCSF;
    }
    _traversalTimes[_currentConfig] = time;
  }

  inline long getEvidence(Configuration configuration) const override { return _traversalTimes.at(configuration); }

  inline void reset(size_t iteration) override {
    _invalidConfigs.clear();
    _searchSpaceSizeNoCSF = 0;
    _numTestedConfigsNoCSF = 0;
    _traversalTimes.clear();
    tune();
  }

  inline bool tune(bool currentInvalid = false) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override { return _containerOptions; }

  inline bool searchSpaceIsTrivial() const override;

  inline bool searchSpaceIsEmpty() const override;

 private:
  inline void selectOptimalConfiguration();

  std::set<ContainerOption> _containerOptions;
  std::set<TraversalOption> _traversalOptions;
  std::set<LoadEstimatorOption> _loadEstimatorOptions;
  std::set<DataLayoutOption> _dataLayoutOptions;
  std::set<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;

  Configuration _currentConfig;
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
  std::set<Configuration> _invalidConfigs;

  Random _rng;
  size_t _maxEvidence;
  size_t _searchSpaceSizeNoCSF;
  size_t _numTestedConfigsNoCSF;
};

bool RandomSearch::tune(bool currentInvalid) {
  // We can overwrite the cellSizeFactor, because _currentConfig will be changed anyways.
  if (currentInvalid) {
    _currentConfig.cellSizeFactor = 1;
    if (_invalidConfigs.find(_currentConfig) == _invalidConfigs.end()) {
      ++_numTestedConfigsNoCSF;
    }
    _invalidConfigs.emplace(_currentConfig);
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = Configuration();
    return false;
  }

  if (_traversalTimes.size() >= _maxEvidence || _numTestedConfigsNoCSF >= _searchSpaceSizeNoCSF) {
    // select best config
    selectOptimalConfiguration();
    return false;
  }

  // needed to test if a new configuration has been found for infinite cellSizeFactors
  Configuration testingConfig;
  do {
    // select random config
    _currentConfig.container = _rng.pickRandom(_containerOptions);
    _currentConfig.cellSizeFactor = _cellSizeFactors->getRandom(_rng);
    _currentConfig.traversal = _rng.pickRandom(_traversalOptions);
    auto applicableLoadEstimators = loadEstimators::getApplicableLoadEstimators(
        _currentConfig.container, _currentConfig.traversal, _loadEstimatorOptions);
    _currentConfig.loadEstimator = _rng.pickRandom(applicableLoadEstimators);
    _currentConfig.dataLayout = _rng.pickRandom(_dataLayoutOptions);
    _currentConfig.newton3 = _rng.pickRandom(_newton3Options);

    testingConfig = Configuration(_currentConfig);
    if (not _cellSizeFactors->isFinite()) {
      testingConfig.cellSizeFactor = 1;
    }
  } while (_invalidConfigs.find(testingConfig) != _invalidConfigs.end() ||
           _traversalTimes.find(_currentConfig) != _traversalTimes.end());
  return true;
}

void RandomSearch::selectOptimalConfiguration() {
  if (searchSpaceIsTrivial()) {
    _currentConfig.container = *_containerOptions.begin();
    _currentConfig.cellSizeFactor = _cellSizeFactors->getMin();
    _currentConfig.traversal = *_traversalOptions.begin();
    auto applicableLoadEstimators = loadEstimators::getApplicableLoadEstimators(
        _currentConfig.container, _currentConfig.traversal, _loadEstimatorOptions);
    _currentConfig.loadEstimator = *applicableLoadEstimators.begin();
    _currentConfig.dataLayout = *_dataLayoutOptions.begin();
    _currentConfig.newton3 = *_newton3Options.begin();
    AutoPasLog(debug, "Selected optimal configuration {}", _currentConfig.toString());
    return;
  }

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "RandomSearch: Trying to determine fastest configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  auto optimum = std::min_element(_traversalTimes.begin(), _traversalTimes.end(),
                                  [](std::pair<Configuration, size_t> a, std::pair<Configuration, size_t> b) -> bool {
                                    return a.second < b.second;
                                  });

  _currentConfig = optimum->first;

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig.toString());
}

bool RandomSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  // if no load estimators are specified, none is automatically added.
  return _containerOptions.size() == 1 and (_cellSizeFactors->isFinite() && _cellSizeFactors->size() == 1) and
         _traversalOptions.size() == 1 and _loadEstimatorOptions.size() <= 1 and _dataLayoutOptions.size() == 1 and
         _newton3Options.size() == 1;
}

bool RandomSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerOptions.empty() or (_cellSizeFactors->isFinite() && _cellSizeFactors->size() == 0) or
         _traversalOptions.empty() or _dataLayoutOptions.empty() or _newton3Options.empty();
}

void RandomSearch::removeN3Option(Newton3Option badNewton3Option) {
  _newton3Options.erase(badNewton3Option);

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

}  // namespace autopas
