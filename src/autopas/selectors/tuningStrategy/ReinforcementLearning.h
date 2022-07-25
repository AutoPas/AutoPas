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
#include <cstdlib>

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

    double _oldState = 0;

    auto _Entry = _states.find(*_currentConfig);
    //        maybe change to pre initialize states to reduce places were error can occur
    if (_Entry != _states.end()) {
      _oldState = _Entry->second;
    }

    //    if(iteration>0){
    //      _oldState = _states.at(*_currentConfig);
    //    }

    //    gamma determines the importance of future rewards
    //    without gamma : gamma = 1
    //    sets a lot of importance on the first run/initializing
    //    double _newState = _oldState + alpha * (time-_oldState);

    //    with gamma
    //    double _newState = _oldState + _alpha * (- time + (_gamma * _oldState - _oldState));

    double _newState = 0;
    if (_firstTuningPhase) {
      _newState = -time;
    } else {

//      gamma
      _newState = _oldState + _alpha * (-time - _gamma * _oldState);
//        no gamma
//      _newState = _oldState + _alpha * (-time - _oldState);
    }

    _states[*_currentConfig] = _newState;
  }

  inline long getEvidence(Configuration configuration) const override { return _traversalTimes.at(configuration); }

  inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  inline void reset(size_t iteration) override {
    _traversalTimes.clear();
    if (_firstTuningPhase) {
      _currentConfig = _searchSpace.begin();
    } else {
      _selectedConfigs = _getCollectionOfConfigurations();
      _currentConfig = _selectedConfigs.begin();
    }
  }

  inline bool tune(bool = false) override;

  inline double getState(Configuration configuration) const { return _states.at(configuration); }

  inline bool getFirstTuningPhase() const { return _firstTuningPhase; }

  inline bool getStartTuningPhase() const { return _startTuningPhase; }

  inline std::set<Configuration> getSearchSpace() const { return _searchSpace; }

  std::unordered_map<Configuration, double, ConfigHash> _states;

 private:
  inline void selectOptimalConfiguration();

  inline std::set<Configuration> _getCollectionOfConfigurations();

  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;

  bool _firstTuningPhase = true;
  bool _startTuningPhase = true;
  const char * _alpha_tmp  = std::getenv("ALPHA");
  double _alpha = (_alpha_tmp== nullptr) ? 1 : atof(_alpha_tmp);
  const char * _gamma_tmp = std::getenv("GAMMA");
  double _gamma = (_gamma_tmp==nullptr) ? 1 : atof(_gamma_tmp);
  size_t _randomExplorations = 5;
  std::set<Configuration> _selectedConfigs{};
};

bool ReinforcementLearning::tune(bool) {
  // repeat as long as traversals are not applicable or we run out of configs
  if (_firstTuningPhase) {
    // run everything once
    ++_currentConfig;
    if (_currentConfig == _searchSpace.end()) {
      selectOptimalConfiguration();
      _firstTuningPhase = false;
      return false;
    }

    return true;
  } else {
    if (_startTuningPhase) {
      _startTuningPhase = false;
    }

    ++_currentConfig;

    if (_currentConfig == _selectedConfigs.end()) {
      selectOptimalConfiguration();
      _startTuningPhase = true;
      return false;
    }
    return true;
  }
}

std::set<Configuration> ReinforcementLearning::_getCollectionOfConfigurations() {
  if (_randomExplorations >= _searchSpace.size()) {
    return _searchSpace;
  }

  std::set<Configuration> _collection;
  _collection.insert(*_currentConfig);

  int _random = _randomExplorations;
  while (_random > 0) {
    size_t _size = _searchSpace.size() - 1;

    std::random_device rd;                            // obtain a random number from hardware
    std::mt19937 gen(rd());                           // seed the generator
    std::uniform_int_distribution<> distr(0, _size);  // define the range
    size_t _place = distr(gen);                       // generate numbers

    if (_place >= _searchSpace.size()) {
      autopas::utils::ExceptionHandler::exception("Reinforcement Learning: Selected Number too high!");
    }
    auto _selectedConfig = _searchSpace.begin();
    std::advance(_selectedConfig, _place);
    auto res = _collection.insert(*_selectedConfig);
    if (res.second) {
      _random--;
    }
  }
  if (_collection.size() <= 0) {
    autopas::utils::ExceptionHandler::exception("Reinforcement Learning: Collection is empty!");
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

  if (_searchSpace.empty()) {
    utils::ExceptionHandler::exception("Reinforcement Learning: Search Space is empty");
  }

  const auto optimum = std::max_element(_states.begin(), _states.end(),
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
