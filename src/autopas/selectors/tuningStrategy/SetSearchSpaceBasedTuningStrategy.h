/**
 * @file SetSearchSpaceBasedTuningStrategy.h
 * @author Julian Pelloth
 * @date 29.04.2020
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Super class for tuning strategies with a set based search space.
 */
class SetSearchSpaceBasedTuningStrategy : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the SetSearchSpaceBasedTuningStrategy that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param allowedVerletRebuildFrequencies
   */
  SetSearchSpaceBasedTuningStrategy(const std::set<ContainerOption> &allowedContainerOptions,
                                    const std::set<double> &allowedCellSizeFactors,
                                    const std::set<TraversalOption> &allowedTraversalOptions,
                                    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                    const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                    const std::set<Newton3Option> &allowedNewton3Options,
                                    const std::set<int> &allowedVerletRebuildFrequencies)
      : _allowedContainerOptions(allowedContainerOptions) {
    // sets search space and current config
    populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                        allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options, allowedVerletRebuildFrequencies);
  }

  /**
   * Constructor for the SetSearchSpaceBasedTuningStrategy that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit SetSearchSpaceBasedTuningStrategy(std::set<Configuration> allowedConfigurations)
      : _allowedContainerOptions{}, _searchSpace(std::move(allowedConfigurations)) {
    for (const auto &config : _searchSpace) {
      _allowedContainerOptions.insert(config.container);
    }
    _currentConfig = _searchSpace.begin();
  }

  [[nodiscard]] inline std::set<ContainerOption> getAllowedContainerOptions() const override {
    return _allowedContainerOptions;
  }

  inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  [[nodiscard]] inline bool searchSpaceIsTrivial() const override { return _searchSpace.size() == 1; }

  [[nodiscard]] inline bool searchSpaceIsEmpty() const override { return _searchSpace.empty(); }

  [[nodiscard]] inline bool smoothedHomogeneityAndMaxDensityNeeded() const override { return false; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

 protected:
  /**
   * Fills the search space with the cartesian product of the given options (minus invalid combinations).
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedVerletRebuildFrequencies
   */
  inline void populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                  const std::set<double> &allowedCellSizeFactors,
                                  const std::set<TraversalOption> &allowedTraversalOptions,
                                  const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options,
                                  const std::set<int> &allowedVerletRebuildFrequencies);

  /**
   * Finds the optimal configuration in a given search space.
   * @param searchSpace Map<Configuration, runTime> to search.
   * @return Iterator to optimal (=fastest) configuration.
   */
  static inline auto getOptimum(const std::unordered_map<Configuration, long, ConfigHash> &searchSpace);

  /**
   * Contains every allowed container option.
   */
  std::set<ContainerOption> _allowedContainerOptions{};
  /**
   * Contains every configuration.
   */
  std::set<Configuration> _searchSpace{};
  /**
   * Contains the corrent configuration.
   */
  std::set<Configuration>::iterator _currentConfig{};
};

void SetSearchSpaceBasedTuningStrategy::populateSearchSpace(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
    const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options,
    const std::set<int> &allowedVerletRebuildFrequencies) {
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
        // if load estimators are not applicable LoadEstimatorOption::none is returned.
        const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
            loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
        for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
          for (const auto &dataLayoutOption : allowedDataLayoutOptions) {
            for (const auto &newton3Option : allowedNewton3Options) {
              for (const auto &verletRebuildFrequency : allowedVerletRebuildFrequencies) {
                _searchSpace.emplace(containerOption, cellSizeFactor, traversalOption, loadEstimatorOption,
                                     dataLayoutOption, newton3Option, verletRebuildFrequency);
              }
            }
          }
        }
      }
  }

  AutoPasLog(debug, "Points in search space: {}", _searchSpace.size());

  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("No valid configurations could be created.");
  }

  _currentConfig = _searchSpace.begin();
}

auto SetSearchSpaceBasedTuningStrategy::getOptimum(
    const std::unordered_map<Configuration, long, ConfigHash> &searchSpace) {
  // find mapping with smallest second (=runtime)
  return std::min_element(
      searchSpace.begin(), searchSpace.end(),
      [](std::pair<Configuration, long> a, std::pair<Configuration, long> b) -> bool { return a.second < b.second; });
}

void SetSearchSpaceBasedTuningStrategy::removeN3Option(Newton3Option badNewton3Option) {
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
