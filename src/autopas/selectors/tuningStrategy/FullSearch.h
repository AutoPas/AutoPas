/**
 * @file FullSearch.h
 * @author F. Gratl
 * @date 5/29/19
 */

#pragma once

#include <set>
#include <sstream>
#include "TuningStrategyInterface.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

class FullSearch : public TuningStrategyInterface {
 public:
  FullSearch(const std::set<ContainerOption> &allowedContainerOptions,
             const std::set<TraversalOption> &allowedTraversalOptions,
             const std::set<DataLayoutOption> &allowedDataLayoutOptions,
             const std::set<Newton3Option> &allowedNewton3Options, SelectorStrategy selectorStrategy)
      : _selectorStrategy(selectorStrategy), _containerOptions(allowedContainerOptions) {
    // sets search space and current config
    generateSearchSpace(allowedContainerOptions, allowedTraversalOptions, allowedDataLayoutOptions,
                        allowedNewton3Options);
  }

  /**
   * Constructor for the FullSearch that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   * @param selectorStrategy Strategy for the configuration selection.
   */
  FullSearch(const std::set<Configuration> &allowedConfigurations, SelectorStrategy selectorStrategy)
      : _selectorStrategy(selectorStrategy),
        _containerOptions(),
        _searchSpace(allowedConfigurations),
        _currentConfig(_searchSpace.begin()) {
    for (auto config : _searchSpace) {
      _containerOptions.insert(config._container);
    }
  }

  inline Configuration getCurrentConfiguration() override { return *_currentConfig; }

  inline void addEvidence(long time) override { _traversalTimes[*_currentConfig].push_back(time); }

  inline void reset() override {
    _traversalTimes.clear();
    _currentConfig = _searchSpace.begin();
  }

  inline bool tune() override;

  inline std::set<ContainerOption> getAllowedContainerOptions() override { return _containerOptions; };

  inline bool searchSpaceOneOption() override { return _searchSpace.size() == 1; }

  inline bool searchSpaceEmpty() override { return _searchSpace.empty(); }

 private:
  inline void generateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                  const std::set<TraversalOption> &allowedTraversalOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options);

  inline void selectOptimalConfiguration();

  SelectorStrategy _selectorStrategy;
  std::set<ContainerOption> _containerOptions;
  std::set<Configuration> _searchSpace;
  std::set<Configuration>::iterator _currentConfig;
  std::unordered_map<Configuration, std::vector<size_t>, ConfigHash> _traversalTimes;
};

void FullSearch::generateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                     const std::set<TraversalOption> &___allowedTraversalOptions,
                                     const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                     const std::set<Newton3Option> &allowedNewton3Options) {
  auto dummySet = {TraversalOption::dummyTraversal};
  std::set<TraversalOption> allowedTraversalOptions;
  std::set_union(___allowedTraversalOptions.begin(), ___allowedTraversalOptions.end(), dummySet.begin(), dummySet.end(),
                 std::inserter(allowedTraversalOptions,allowedTraversalOptions.begin()));

  // generate all potential configs
  for (auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    std::set<TraversalOption> allContainerTraversals = applicableTraversals::allApplicableTraversals(containerOption);
    //@TODO dummyTraversal needed until all containers support propper traversals
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (auto &traversalOption : allowedAndApplicable) {
      for (auto &dataLayoutOption : allowedDataLayoutOptions) {
        for (auto &newton3Option : allowedNewton3Options) {
          _searchSpace.emplace(containerOption, traversalOption, dataLayoutOption, newton3Option);
        }
      }
    }
  }

  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("FullSearch: No valid configurations could be created.");
  }

  _currentConfig = _searchSpace.begin();
}

bool FullSearch::tune() {
  // repeat as long as traversals are not applicable or we run out of configs
  if (_currentConfig != _searchSpace.end()) {
    ++_currentConfig;
  } else {  // reached end of tuning phase
    // sets _currentConfig
    selectOptimalConfiguration();
    return false;
  }

  return true;
}

void FullSearch::selectOptimalConfiguration() {
  if (_searchSpace.size() == 1) {
    _currentConfig = _searchSpace.begin();
    return;
  }

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "FullSearch: Trying to determine fastest configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  long optimalTraversalTime = std::numeric_limits<long>::max();
  Configuration optimalConfiguration;
  // reduce sample values
  for (auto &configAndTimes : _traversalTimes) {
    long value = OptimumSelector::optimumValue(configAndTimes.second, _selectorStrategy);

    // save all values for debugging purposes
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      configAndTimes.second.push_back(value);
    }

    if (value < optimalTraversalTime) {
      optimalTraversalTime = value;
      optimalConfiguration = configAndTimes.first;
    }
  }

  // print all configs, times and their reduced values
  if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
    std::stringstream ss;
    // print all configs
    for (auto &p : _traversalTimes) {
      ss << std::endl << p.first.toString() << " : [";
      // print all timings
      for (size_t i = 0; i < p.second.size() - 1; ++i) {
        ss << " " << p.second[i];
      }
      ss << " ] ";
      ss << "Reduced value: " << p.second[p.second.size() - 1];
    }
    AutoPasLog(debug, "Collected times: {}", ss.str());
  }

  _currentConfig = _searchSpace.find(optimalConfiguration);
  // sanity check
  if (_currentConfig == _searchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "FullSearch: Optimal configuration not found in list of configurations!");
  }

  // measurements are not needed anymore
  _traversalTimes.clear();

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}

}  // namespace autopas
