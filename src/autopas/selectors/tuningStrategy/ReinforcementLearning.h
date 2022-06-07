/**
 * @file ReinforcementLearning.h
 * @author L. Laumeyer
 * @date 17.05.2022
 */

#pragma once

#include <random>
#include <set>
#include <sstream>
#include <utility>

#include "SetSearchSpaceBasedTuningStrategy.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Using reinforcement learning to determine the best configuration
 */
class ReinforcementLearning : public SetSearchSpaceBasedTuningStrategy {
 public:
  /**
   * Constructor for the Reinforcement Learning Tuning Strategy that generates the search space from the allowed
   * options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   */
  ReinforcementLearning(const std::set<ContainerOption> &allowedContainerOptions,
                        const std::set<double> &allowedCellSizeFactors,
                        const std::set<TraversalOption> &allowedTraversalOptions,
                        const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                        const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                        const std::set<Newton3Option> &allowedNewton3Options)
      : SetSearchSpaceBasedTuningStrategy(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                                          allowedLoadEstimatorOptions, allowedDataLayoutOptions,
                                          allowedNewton3Options) {}

  /**
   * Constructor for the ReinforcementLearning that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit ReinforcementLearning(std::set<Configuration> allowedConfigurations)
      : SetSearchSpaceBasedTuningStrategy(std::move(allowedConfigurations)) {}

  inline void addEvidence(long time, size_t iteration) override {
    //  use time as a reward would require to look for smallest reward
    _traversalTimes[*_currentConfig] = time;
    auto _Entry = _traversalTimes.find(*_currentConfig);
    double _oldState = 0;
    //    maybe change to pre initialize states to reduce places were error can occure
    if (_Entry != _traversalTimes.end()) {
      _oldState = _Entry->second;
    }
    //    gamma determines the importance of future rewards
    //    without gamma : gamma = 1
    //    sets a lot of importance on the first run/initializing
    //    double _newState = _oldState + alpha * (time-_oldState);

    //    with gamma
    double _newState = _oldState + _alpha * (time + (_gamma * _oldState - _oldState));

    _states[*_currentConfig] = _newState;
  }

  inline long getEvidence(Configuration configuration) const override { return _traversalTimes.at(configuration); }

  inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  inline void reset(size_t iteration) override { _traversalTimes.clear(); }

  inline bool tune(bool = false) override;

 private:
  inline void selectOptimalConfiguration();

  inline std::set<Configuration> _getCollectionOfConfigurations();

  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
  std::unordered_map<Configuration, double, ConfigHash> _states;
  bool _firstRun = true;
  double _alpha = 0.9;
  double _gamma = 1;
  size_t _randomExplorations = 5;
};

bool ReinforcementLearning::tune(bool) {
  // repeat as long as traversals are not applicable or we run out of configs
  if (_firstRun) {
    // run everything once
    ++_currentConfig;
    if (_currentConfig == _searchSpace.end()) {
      selectOptimalConfiguration();
      _firstRun = false;
      return false;
    }

    return true;
  } else {
    std::set<Configuration> _selectedConfigs{};
    _selectedConfigs = _getCollectionOfConfigurations();
    _currentConfig = _selectedConfigs.begin();
    ++_currentConfig;

    if (_currentConfig == _selectedConfigs.end()) {
      selectOptimalConfiguration();
      _firstRun = false;
      return false;
    }
    return true;
  }
}

std::set<Configuration> ReinforcementLearning::_getCollectionOfConfigurations() {
  std::set<Configuration> _collection;
  selectOptimalConfiguration();
  _collection.insert(*_currentConfig);

  if (_randomExplorations >= _searchSpace.size()) {
    return _searchSpace;
  }
  int _random = _randomExplorations;
  while (_random > 0) {
    size_t _size = _searchSpace.size();

    std::random_device rd;                            // obtain a random number from hardware
    std::mt19937 gen(rd());                           // seed the generator
    std::uniform_int_distribution<> distr(0, _size);  // define the range
    int _place = distr(gen);                          // generate numbers

    auto _selectedConfig = _searchSpace.begin();
    std::advance(_selectedConfig, _place);
    auto res = _collection.insert(*_selectedConfig);
    if (res.second) {
      _random--;
    }
  }
  return _collection;
}

void ReinforcementLearning::selectOptimalConfiguration() {
  if (_searchSpace.size() == 1) {
    _currentConfig = _searchSpace.begin();
    return;
  }

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "Reinforcement Learning: Trying to determine the best configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  const auto optimum = std::min_element(_states.begin(), _states.end(),
                                        [](std::pair<Configuration, double> a,
                                           std::pair<Configuration, double> b) -> bool { return a.second < b.second; });

  _currentConfig = _searchSpace.find(optimum->first);
  // sanity check
  if (_currentConfig == _searchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "Reinforcement Learning: Optimal configuration not found in list of configurations!");
  }
}

}  // namespace autopas
