/**
 * @file TuningStrategySuperClass.h
 * @author Julian Pelloth
 * @date 29.04.2020
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Super class for tuning strategies
 */
class TuningStrategySuperClass : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the TuningStrategySuperClass that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   */
  TuningStrategySuperClass(const std::set<ContainerOption> &allowedContainerOptions,
                           const std::set<double> &allowedCellSizeFactors,
                           const std::set<TraversalOption> &allowedTraversalOptions,
                           const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                           const std::set<Newton3Option> &allowedNewton3Options)
      : _containerOptions(allowedContainerOptions) {
    // sets search space and current config
    populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                        allowedDataLayoutOptions, allowedNewton3Options);
  }

  /**
   * Constructor for the TuningStrategySuperClass that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit TuningStrategySuperClass(std::set<Configuration> allowedConfigurations)
      : _containerOptions{}, _searchSpace(std::move(allowedConfigurations)) {
    for (const auto &config : _searchSpace) {
      _containerOptions.insert(config.container);
    }
  }

  inline std::set<ContainerOption> getAllowedContainerOptions() const override { return _containerOptions; }

  inline bool searchSpaceIsTrivial() const override { return _searchSpace.size() == 1; }

  inline bool searchSpaceIsEmpty() const override { return _searchSpace.empty(); }

 protected:
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

  /**
   * Finds the optimum in an unordered map
   * @param findOptimumSearchSpace
   * @return optimum
   */
  static inline std::__detail::_Node_const_iterator<std::pair<const Configuration, unsigned long>, false, true>
  getOptimum(const std::unordered_map<Configuration, size_t, ConfigHash> &findOptimumSearchSpace);

  std::set<ContainerOption> _containerOptions;
  std::set<Configuration> _searchSpace;
};

void TuningStrategySuperClass::populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
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

std::__detail::_Node_const_iterator<std::pair<const Configuration, unsigned long>, false, true>
TuningStrategySuperClass::getOptimum(
    const std::unordered_map<Configuration, size_t, ConfigHash> &findOptimumSearchSpace) {
  return std::min_element(findOptimumSearchSpace.begin(), findOptimumSearchSpace.end(),
                          [](std::pair<Configuration, size_t> a, std::pair<Configuration, size_t> b) -> bool {
                            return a.second < b.second;
                          });
}
}  // namespace autopas