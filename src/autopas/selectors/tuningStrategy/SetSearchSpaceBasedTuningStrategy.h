/**
 * @file SetSearchSpaceBasedTuningStrategy.h
 * @author Julian Pelloth
 * @date 29.04.2020
 */

#pragma once

#include "TuningStrategyInterface.h"
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
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   */
  SetSearchSpaceBasedTuningStrategy(const std::set<ContainerOption> &allowedContainerOptions,
                                    const std::set<double> &allowedCellSizeFactors,
                                    const std::set<TraversalOption> &allowedTraversalOptions,
                                    const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                    const std::set<Newton3Option> &allowedNewton3Options)
      : _allowedContainerOptions(allowedContainerOptions) {
    // sets search space and current config
    populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                        allowedDataLayoutOptions, allowedNewton3Options);
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
  }

  [[nodiscard]] inline std::set<ContainerOption> getAllowedContainerOptions() const override {
    return _allowedContainerOptions;
  }

  [[nodiscard]] inline bool searchSpaceIsTrivial() const override { return _searchSpace.size() == 1; }

  [[nodiscard]] inline bool searchSpaceIsEmpty() const override { return _searchSpace.empty(); }

 protected:
  /**
   * Fills the search space with the cartesian product of the given options (minus invalid combinations).
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  inline void populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                  const std::set<double> &allowedCellSizeFactors,
                                  const std::set<TraversalOption> &allowedTraversalOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options);

  /**
   * Finds the optimum in an unordered map
   * @param findOptimumSearchSpace
   * @return optimum
   */
  static inline auto getOptimum(const std::unordered_map<Configuration, size_t, ConfigHash> &findOptimumSearchSpace);

  /**
   * Contains every allowed container option.
   */
  std::set<ContainerOption> _allowedContainerOptions;
  /**
   * Contains every configuration.
   */
  std::set<Configuration> _searchSpace;
};

void SetSearchSpaceBasedTuningStrategy::populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                                            const std::set<double> &allowedCellSizeFactors,
                                                            const std::set<TraversalOption> &allowedTraversalOptions,
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
        for (const auto &dataLayoutOption : allowedDataLayoutOptions) {
          for (const auto &newton3Option : allowedNewton3Options) {
            _searchSpace.emplace(containerOption, cellSizeFactor, traversalOption, dataLayoutOption, newton3Option);
          }
        }
      }
  }

  AutoPasLog(debug, "Points in search space: {}", _searchSpace.size());

  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("FullSearch: No valid configurations could be created.");
  }
}

auto SetSearchSpaceBasedTuningStrategy::getOptimum(
    const std::unordered_map<Configuration, size_t, ConfigHash> &findOptimumSearchSpace) {
  return std::min_element(findOptimumSearchSpace.begin(), findOptimumSearchSpace.end(),
                          [](std::pair<Configuration, size_t> a, std::pair<Configuration, size_t> b) -> bool {
                            return a.second < b.second;
                          });
}
}  // namespace autopas