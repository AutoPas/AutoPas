/**
 * @file GenerateValidConfigurations.h
 * @author The AI ghost of F. Gratl
 * @date 24.03.26
 */

#pragma once

#include <set>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/NumberSetFinite.h"

/**
 * Struct to hold a container and a cell size factor.
 */
struct ContainerConfiguration {
  autopas::ContainerOption container;
  double cellSizeFactor;

  /**
   * Equality operator.
   * @param rhs
   * @return
   */
  bool operator==(const ContainerConfiguration &rhs) const {
    return container == rhs.container and cellSizeFactor == rhs.cellSizeFactor;
  }

  /**
   * Inequality operator.
   * @param rhs
   * @return
   */
  bool operator!=(const ContainerConfiguration &rhs) const { return not(*this == rhs); }

  /**
   * Comparison operator.
   * @param rhs
   * @return
   */
  bool operator<(const ContainerConfiguration &rhs) const {
    return std::tie(container, cellSizeFactor) < std::tie(rhs.container, rhs.cellSizeFactor);
  }

  /**
   * String representation of the configuration.
   * @return
   */
  [[nodiscard]] std::string toString() const {
    return "{Container: " + container.to_string() + " , CellSizeFactor: " + std::to_string(cellSizeFactor) + "}";
  }
};

/**
 * Stream operator for ContainerConfiguration.
 * @param os
 * @param configuration
 * @return
 */
inline std::ostream &operator<<(std::ostream &os, const ContainerConfiguration &configuration) {
  return os << configuration.toString();
}

/**
 * Generates all valid configurations for a given interaction type or all interaction types.
 *
 * Intended to be used with the default arguments. If one is overridden, make sure you have a good reason and write down
 * why.
 *
 * @param interactionType The interaction type. If provided with "all", all valid configurations for all interaction
 * types are returned.
 * @param allowedContainerOptions By default, all options.
 * @param allowedTraversalOptions By default, all options.
 * @param allowedLoadEstimatorOptions By default, all options.
 * @param allowedDataLayoutOptions By default, all options.
 * @param allowedNewton3Options By default, all options.
 * @param allowedCellSizeFactors By default, {0.5, 1.0, 1.5}
 * @return
 */
inline std::set<autopas::Configuration> generateAllValidConfigurations(
    autopas::InteractionTypeOption interactionType,
    const std::set<autopas::ContainerOption> &allowedContainerOptions = autopas::ContainerOption::getAllOptions(),
    const std::set<autopas::TraversalOption> &allowedTraversalOptions = autopas::TraversalOption::getAllOptions(),
    const std::set<autopas::LoadEstimatorOption> &allowedLoadEstimatorOptions =
        autopas::LoadEstimatorOption::getAllOptions(),
    const std::set<autopas::DataLayoutOption> &allowedDataLayoutOptions = autopas::DataLayoutOption::getAllOptions(),
    const std::set<autopas::Newton3Option> &allowedNewton3Options = autopas::Newton3Option::getAllOptions(),
    const std::set<double> &allowedCellSizeFactors = {0.5, 1.0, 1.5}) {
  if (interactionType == autopas::InteractionTypeOption::all) {
    std::set<autopas::Configuration> allConfigs;
    for (auto iType : autopas::InteractionTypeOption::getMostOptions()) {
      const autopas::NumberSetFinite<double> csfs(allowedCellSizeFactors);
      const auto configs = autopas::SearchSpaceGenerators::cartesianProduct(
          allowedContainerOptions, allowedTraversalOptions, allowedLoadEstimatorOptions, allowedDataLayoutOptions,
          allowedNewton3Options, &csfs, iType);
      allConfigs.insert(configs.begin(), configs.end());
    }
    return allConfigs;
  } else {
    const autopas::NumberSetFinite<double> csfs(allowedCellSizeFactors);
    return autopas::SearchSpaceGenerators::cartesianProduct(allowedContainerOptions, allowedTraversalOptions,
                                                            allowedLoadEstimatorOptions, allowedDataLayoutOptions,
                                                            allowedNewton3Options, &csfs, interactionType);
  }
}

/**
 * Generates all valid container configurations.
 * Intended to be used with the default arguments. If one is overridden, make sure you have a good reason and write down
 * why.
 * @param allowedContainerOptions By default, all options.
 * @param allowedCellSizeFactors By default, {0.5, 1.0, 1.5}
 * @return
 */
inline std::set<ContainerConfiguration> generateAllValidContainerConfigurations(
    const std::set<autopas::ContainerOption> &allowedContainerOptions = autopas::ContainerOption::getAllOptions(),
    const std::set<double> &allowedCellSizeFactors = {0.5, 1.0, 1.5}) {
  std::set<ContainerConfiguration> containerConfigs;
  for (const auto &containerOption : allowedContainerOptions) {
    for (const auto csf : allowedCellSizeFactors) {
      // Create a dummy configuration to check validity
      // We use pairwise interaction as default, since it should not matter for container/csf compatibility.
      // We also use the first valid traversal for this container.
      const auto interactionType = autopas::InteractionTypeOption::pairwise;
      const auto traversals = autopas::compatibleTraversals::allCompatibleTraversals(containerOption, interactionType);
      if (traversals.empty()) {
        continue;
      }
      const auto &traversalOption = *traversals.begin();
      const auto loadEstimators = autopas::loadEstimators::getApplicableLoadEstimators(
          containerOption, traversalOption, autopas::LoadEstimatorOption::getAllOptions());
      const auto &loadEstimatorOption = *loadEstimators.begin();

      const autopas::Configuration configuration{containerOption,
                                                 csf,
                                                 traversalOption,
                                                 loadEstimatorOption,
                                                 autopas::DataLayoutOption::aos,
                                                 autopas::Newton3Option::enabled,
                                                 interactionType};
      if (configuration.hasCompatibleValues()) {
        containerConfigs.insert({containerOption, csf});
      }
    }
  }
  return containerConfigs;
}
