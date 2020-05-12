/**
 * @file FullSearch.h
 * @author F. Gratl
 * @date 29.05.2019
 */

#pragma once

#include <set>
#include <sstream>

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/loadEstimators.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Exhaustive full search of the search space by testing every applicable configuration and then selecting the optimum.
 */
class FullSearch : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the FullSearch that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   */
  FullSearch(const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
             const std::set<TraversalOption> &allowedTraversalOptions,
             const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
             const std::set<DataLayoutOption> &allowedDataLayoutOptions,
             const std::set<Newton3Option> &allowedNewton3Options)
      : _containerOptions(allowedContainerOptions) {
    // sets search space and current config
    populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                        allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options);
  }

  /**
   * Constructor for the FullSearch that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit FullSearch(std::set<Configuration> allowedConfigurations)
      : _containerOptions{}, _searchSpace(std::move(allowedConfigurations)), _currentConfig(_searchSpace.begin()) {
    for (const auto config : _searchSpace) {
      _containerOptions.insert(config.container);
    }
  }

  inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t iteration) override { _traversalTimes[*_currentConfig] = time; }

  inline void reset() override {
    _traversalTimes.clear();
    _currentConfig = _searchSpace.begin();
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
                                  const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options);

  inline void selectOptimalConfiguration();

  std::set<ContainerOption> _containerOptions;
  std::set<Configuration> _searchSpace;
  std::set<Configuration>::iterator _currentConfig;
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
};

void FullSearch::populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                     const std::set<double> &allowedCellSizeFactors,
                                     const std::set<TraversalOption> &allowedTraversalOptions,
                                     const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                     const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                     const std::set<Newton3Option> &allowedNewton3Options) {
  // generate all potential configs
  for (const auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (const auto &cellSizeFactor : allowedCellSizeFactors)
      for (const auto &traversalOption : allowedAndApplicable) {
        const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
            loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
        for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
          for (const auto &dataLayoutOption : allowedDataLayoutOptions) {
            for (const auto &newton3Option : allowedNewton3Options) {
              _searchSpace.emplace(containerOption, cellSizeFactor, traversalOption, loadEstimatorOption,
                                   dataLayoutOption, newton3Option);
            }
          }
        }
      }

    AutoPasLog(debug, "Points in search space: {}", _searchSpace.size());

    if (_searchSpace.empty()) {
      autopas::utils::ExceptionHandler::exception("FullSearch: No valid configurations could be created.");
    }

    _currentConfig = _searchSpace.begin();
  }
}

bool FullSearch::tune(bool) {
  // repeat as long as traversals are not applicable or we run out of configs
  ++_currentConfig;
  if (_currentConfig == _searchSpace.end()) {
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

  const auto optimum = std::min_element(_traversalTimes.begin(), _traversalTimes.end(),
                                        [](std::pair<Configuration, size_t> a,
                                           std::pair<Configuration, size_t> b) -> bool { return a.second < b.second; });

  _currentConfig = _searchSpace.find(optimum->first);
  // sanity check
  if (_currentConfig == _searchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "FullSearch: Optimal configuration not found in list of configurations!");
  }

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}

void FullSearch::removeN3Option(Newton3Option badNewton3Option) {
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
