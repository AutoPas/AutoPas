/**
 * @file PredictiveTuning.h
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#pragma once

#include <set>

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * On every tuning phase, this strategy makes a runtime prediction for every configuration and then
 * only tests those which are within a certain range of the best prediction. In the end, the configuration
 * that performed best during testing is selected.
 */
class PredictiveTuning : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the PredictiveTuning that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param iterations
   */
  PredictiveTuning(const std::set<ContainerOption> &allowedContainerOptions,
                   const std::set<double> &allowedCellSizeFactors,
                   const std::set<TraversalOption> &allowedTraversalOptions,
                   const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                   const std::set<Newton3Option> &allowedNewton3Options, unsigned int &iterations)
      : _containerOptions(allowedContainerOptions), _iterations(iterations) {
    // sets search space and current config
    populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                        allowedDataLayoutOptions, allowedNewton3Options);
  }

  /**
   * Constructor for the PredictiveTuning that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   * @param iterations
   */
  explicit PredictiveTuning(std::set<Configuration> allowedConfigurations, unsigned int &iterations)
      : _containerOptions{},
        _searchSpace(std::move(allowedConfigurations)),
        _currentConfig(_searchSpace.begin()),
        _iterations(iterations) {
    for (auto config : _searchSpace) {
      _containerOptions.insert(config.container);
    }
  }

  inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t iteration) override {
    _traversalTimesStorage[*_currentConfig].emplace_back(std::make_pair(iteration, time));
    _lastTest[*_currentConfig] = _tuningIterationsCounter;
  }

  inline void reset() override {
    _configurationPredictions.clear();
    _optimalSearchSpace.clear();
    _longNotTestedSearchSpace.clear();
    // Too expensive?
    _validSearchSpace = _searchSpace;
    _validConfigurationFound = false;
    _iterationBeginTuningPhase = _iterations;

    selectOptimalSearchSpace();
  }

  inline bool tune(bool = false) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override { return _containerOptions; }

  inline bool searchSpaceIsTrivial() const override { return _searchSpace.size() == 1; }

  inline bool searchSpaceIsEmpty() const override { return _searchSpace.empty(); }

 private:
  /**
   * Fills the search space with the cartesian product of the given options (minus invalid combinations).
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  inline void populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                  const std::set<double> &allowedCellSizeFactors,
                                  const std::set<TraversalOption> &allowedTraversalOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options);

  inline void selectOptimalConfiguration();

  /**
   * Selects the configurations that are going to be tested.
   */
  inline void selectOptimalSearchSpace();

  /**
   * Provides different extrapolation methods for the prediction of the traversal time.
   */
  inline void calculatePredictions();

  /**
   * Predicts the traversal time by placing a line trough the last two traversal points and calculating the prediction
   * for the current time.
   */
  inline void linePrediction();

  /**
   * Creates a new optimalSearchSpace if every configuration in the previous one was invalid.
   */
  inline void reselectOptimalSearchSpace();

  std::set<ContainerOption> _containerOptions;
  std::set<Configuration> _searchSpace;
  std::set<Configuration>::iterator _currentConfig;

  /**
   * Stores the traversal times for each configuration.
   * @param Configuration
   * @param Vector with pairs of the iteration and the traversal time
   */
  std::unordered_map<Configuration, std::vector<std::pair<int, size_t>>, ConfigHash> _traversalTimesStorage;

  /**
   * Contains the predicted time for each configuration.
   * @param Configuration
   * @param traversal prediction
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _configurationPredictions;

  /**
   * Contains the configuration that are predicted to be optimal and going to be tested.
   */
  std::set<Configuration> _optimalSearchSpace;

  /**
   * Contains the configuration that have not been tested for a period of time and are going to be tested.
   */
  std::set<Configuration> _longNotTestedSearchSpace;

  /**
   * Contains the configurations that are not invalid.
   */
  std::set<Configuration> _validSearchSpace;

  /**
   * Stores the the last tuning phase a configuration got tested
   * @param Configuration
   * @param last tuning phase
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _lastTest;

  /**
   * Gets incremented after every completed tuning phase.
   * Mainly used for the traversal time storage.
   */
  int _tuningIterationsCounter = 0;

  /**
   * Pointer to the _iterations variable in AutoTuner
   */
  unsigned int &_iterations;

  /**
   * Stores the iteration at the beginning of a tuning phase
   */
  unsigned int _iterationBeginTuningPhase = 0;

  /**
   * Indicates if a valid configuration was found in the _optimalSearchSpace.
   */
  bool _validConfigurationFound = false;

  /**
   * Intermediate storage for the optimal configuration
   */
  std::set<Configuration>::iterator _optimalConfig;

  /**
   * Factor of the range of the optimal configurations for the optimalSearchSpace
   */
  static constexpr double _relativeOptimumRange = 1.2;

  /**
   * After not being tested this number of tuningPhases a configuration is being emplaced in _optimalSearchSpace
   */
  static constexpr int _maxTuningIterationsWithoutTest = 5;
};

void PredictiveTuning::populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                           const std::set<double> &allowedCellSizeFactors,
                                           const std::set<TraversalOption> &allowedTraversalOptions,
                                           const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                           const std::set<Newton3Option> &allowedNewton3Options) {
  // generate all potential configs
  for (auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (const auto &cellSizeFactor : allowedCellSizeFactors)
      for (auto &traversalOption : allowedAndApplicable) {
        for (auto &dataLayoutOption : allowedDataLayoutOptions) {
          for (auto &newton3Option : allowedNewton3Options) {
            _searchSpace.emplace(containerOption, cellSizeFactor, traversalOption, dataLayoutOption, newton3Option);
          }
        }
      }
  }

  AutoPasLog(debug, "Points in search space: {}", _searchSpace.size());

  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("PredictiveTuning: No valid configurations could be created.");
  }

  for (auto &configuration : _searchSpace) {
    std::vector<std::pair<int, size_t>> vector;
    _traversalTimesStorage.emplace(configuration, vector);
  }

  _currentConfig = _searchSpace.begin();
}

void PredictiveTuning::selectOptimalSearchSpace() {
  if (_searchSpace.size() == 1 || _tuningIterationsCounter < 2) {
    _currentConfig = _searchSpace.begin();
    return;
  }

  calculatePredictions();

  const auto optimum = std::min_element(_configurationPredictions.begin(), _configurationPredictions.end(),
                                        [](std::pair<Configuration, size_t> a,
                                           std::pair<Configuration, size_t> b) -> bool { return a.second < b.second; });

  _optimalSearchSpace.emplace(optimum->first);

  // selects configurations that are near the optimal prediction or have not been tested in a certain number of
  // iterations
  for (auto &configuration : _searchSpace) {
    // Adds configurations that have not been tested for _maxTuningIterationsWithoutTest or are within the
    // _relativeOptimumRange
    if ((float)_configurationPredictions[configuration] / optimum->second <= _relativeOptimumRange) {
      _optimalSearchSpace.emplace(configuration);
    } else if (_tuningIterationsCounter - _lastTest[configuration] >= _maxTuningIterationsWithoutTest) {
      _longNotTestedSearchSpace.emplace(configuration);
    }
  }

  // sanity check
  if (_optimalSearchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("PredicitveTuning: No possible configuration prediction found!");
  }

  _currentConfig = _optimalSearchSpace.begin();
}

void PredictiveTuning::calculatePredictions() {
  linePrediction();
  /*switch (method) {
  * case line: linePrediction(); break;
  * default:

  }*/
}

void PredictiveTuning::linePrediction() {
  for (auto &configuration : _searchSpace) {
    const auto &vector = _traversalTimesStorage[configuration];
    auto traversal1 = vector[vector.size() - 1];
    auto traversal2 = vector[vector.size() - 2];

    // time1 + (time1 - time2) / (iteration1 - iteration2) / tuningPhase - iteration1)
    _configurationPredictions[configuration] = traversal1.second + (traversal1.second - traversal2.second) /
                                                                       (traversal1.first - traversal2.first) *
                                                                       (_iterations - traversal1.first);
  }
}

void PredictiveTuning::reselectOptimalSearchSpace() {
  // This is placed here, because there are no unnecessary removals
  for (auto &configuration : _optimalSearchSpace) {
    _configurationPredictions.erase(configuration);
    _validSearchSpace.erase(configuration);
  }

  _optimalSearchSpace.clear();

  if (_validSearchSpace.size() == 1) {
    _optimalSearchSpace = _validSearchSpace;
    _currentConfig = _optimalSearchSpace.begin();
    return;
  }

  // checks if there are any configurations left that can be tested
  if (_configurationPredictions.empty()) {
    autopas::utils::ExceptionHandler::exception("Predictive Tuning: No valid configuration could be found");
  }

  const auto optimum = std::min_element(_configurationPredictions.begin(), _configurationPredictions.end(),
                                        [](std::pair<Configuration, size_t> a,
                                           std::pair<Configuration, size_t> b) -> bool { return a.second < b.second; });

  if (_validSearchSpace.count(optimum->first) == 0) {
    autopas::utils::ExceptionHandler::exception("Predictive Tuning: No valid optimal configuration could be found");
  }

  _optimalSearchSpace.emplace(optimum->first);

  // selects configurations that are near the optimal prediction
  for (auto &configuration : _validSearchSpace) {
    // Adds configurations that are within the _relativeOptimumRange
    if ((float)_configurationPredictions[configuration] / optimum->second <= _relativeOptimumRange) {
      _optimalSearchSpace.emplace(configuration);
      _longNotTestedSearchSpace.erase(configuration);
    }
  }

  // sanity check
  if (_optimalSearchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("PredicitveTuning: No possible configuration prediction found!");
  }

  _currentConfig = _optimalSearchSpace.begin();
}

bool PredictiveTuning::tune(bool valid) {
  if (not valid) {
    _validConfigurationFound = true;
  }

  // repeat as long as traversals are not applicable or we run out of configs
  ++_currentConfig;

  if (_currentConfig == _searchSpace.end() || _currentConfig == _optimalSearchSpace.end()) {
    if (_validConfigurationFound) {
      selectOptimalConfiguration();
      if (_longNotTestedSearchSpace.empty()) {
        _currentConfig = _optimalConfig;
        _tuningIterationsCounter++;
        return false;
      } else {
        _currentConfig = _longNotTestedSearchSpace.begin();
        return true;
      }
    } else {
      reselectOptimalSearchSpace();
    }
  } else if (_currentConfig == _longNotTestedSearchSpace.end()) {
    _currentConfig = _optimalConfig;
    _tuningIterationsCounter++;
    return false;
  }

  return true;
}

void PredictiveTuning::selectOptimalConfiguration() {
  if (_optimalSearchSpace.size() == 1) {
    _optimalConfig = _optimalSearchSpace.begin();
    return;
  }
  if (_searchSpace.size() == 1) {
    _optimalConfig = _searchSpace.begin();
    return;
  }

  std::unordered_map<Configuration, size_t, ConfigHash> traversalTimes;
  // for (auto &configuration : _optimalSearchSpace) does not work, but would be way more efficient
  for (auto &configuration : _searchSpace) {
    if (_traversalTimesStorage[configuration].back().first >= _iterationBeginTuningPhase) {
      traversalTimes[configuration] = _traversalTimesStorage[configuration].back().second;
    }
  }

  // Time measure strategy
  if (traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "PredictiveTuning: Trying to determine fastest configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  auto optimum = std::min_element(traversalTimes.begin(), traversalTimes.end(),
                                  [](std::pair<Configuration, size_t> a, std::pair<Configuration, size_t> b) -> bool {
                                    return a.second < b.second;
                                  });

  _optimalConfig = _searchSpace.find(optimum->first);
  // sanity check
  if (_optimalConfig == _searchSpace.end() || _optimalConfig == _optimalSearchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "PredicitveTuning: Optimal configuration not found in list of configurations!");
  }

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}

void PredictiveTuning::removeN3Option(Newton3Option badNewton3Option) {
  for (auto ssIter = _searchSpace.begin(); ssIter != _searchSpace.end();) {
    if (ssIter->newton3 == badNewton3Option) {
      // change current config to the next non-deleted
      if (ssIter == _currentConfig) {
        ssIter = _searchSpace.erase(ssIter);
        _currentConfig = ssIter;
      } else {
        ssIter = _searchSpace.erase(ssIter);
      }
    } else {
      ++ssIter;
    }
  }

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

}  // namespace autopas
