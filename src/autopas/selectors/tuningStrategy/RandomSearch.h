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
        _rng(seed),
        _maxEvidence(maxEvidence),
        _initialMaxEvidence(maxEvidence) {
    if (searchSpaceIsEmpty()) {
      autopas::utils::ExceptionHandler::exception("RandomSearch: No valid configurations could be created.");
    }

    tune();
  }

  inline const Configuration &getCurrentConfiguration() const override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t iteration) override {
    if (_traversalTimes.find(_currentConfig) != _traversalTimes.end()) {
      --_maxEvidence;
    }
    _traversalTimes[_currentConfig] = time;
  }

  inline long getEvidence(Configuration configuration) const override { return _traversalTimes.at(configuration); }

  inline void reset(size_t iteration) override {
    _traversalTimes.clear();
    _maxEvidence = _initialMaxEvidence;
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

  Random _rng;
  size_t _maxEvidence;
  size_t _initialMaxEvidence;
};

bool RandomSearch::tune(bool currentInvalid) {
  if (currentInvalid) {
    --_maxEvidence;
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = Configuration();
    return false;
  }

  if (_traversalTimes.size() >= _maxEvidence) {
    // select best config
    selectOptimalConfiguration();
    return false;
  }

  // select random config
  _currentConfig.container = _rng.pickRandom(_containerOptions);
  _currentConfig.cellSizeFactor = _cellSizeFactors->getRandom(_rng);
  _currentConfig.traversal = _rng.pickRandom(_traversalOptions);
  auto applicableLoadEstimators = loadEstimators::getApplicableLoadEstimators(
      _currentConfig.container, _currentConfig.traversal, _loadEstimatorOptions);
  _currentConfig.loadEstimator = _rng.pickRandom(applicableLoadEstimators);
  _currentConfig.dataLayout = _rng.pickRandom(_dataLayoutOptions);
  _currentConfig.newton3 = _rng.pickRandom(_newton3Options);
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
