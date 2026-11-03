/**
 * @file GenerateValidConfigurations.h
 * @author The AI ghost of F. Gratl
 * @date 24.03.26
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <set>
#include <string>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/OpenMPKindOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/NumberSetFinite.h"

/**
 * OpenMP schedule kinds and chunk sizes to test.
 */
struct OpenMPSchedulesToTest {
  /**
   * The OpenMP schedule kinds.
   */
  std::set<autopas::OpenMPKindOption> kinds;
  /**
   * The OpenMP chunk sizes.
   */
  std::set<autopas::OpenMPChunkSize> chunkSizes;
};

/**
 * Gets the OpenMP schedules that the configuration sweeps test, chosen at runtime with the environment variable
 * AUTOPAS_TEST_OPENMP_SCHEDULES:
 * - static1 (default, also if unset or empty): static with chunk size 1, since OpenMP scheduling is assumed to be
 *   correct and therefore does not need to be swept over. To be sure that this assumption holds, more detailed sweeps
 *   can be obtained with the following env variable options.
 * - core: static, dynamic, and static_steal, with chunk sizes 1, balanced, and golden. The scheduling kinds known to be
 *   potentially optimal with AutoPas, plus the recommended chunk sizes.
 * - extended: static, dynamic, guided, auto, and LB4OMP's scheduling techniques used by Auto4OMP, with chunk sizes 1,
 *   2, 1000 (stress-test), balanced, and golden.
 *
 * Kinds that the OpenMP runtime does not support, and incompatible pairs of kind and chunk size, are filtered out by
 * Configuration::hasCompatibleValues().
 *
 * @return
 */
inline OpenMPSchedulesToTest getOpenMPSchedulesToTest() {
  const char *suiteEnv = std::getenv("AUTOPAS_TEST_OPENMP_SCHEDULES");
  const std::string suite = suiteEnv == nullptr ? "" : suiteEnv;
  if (suite.empty() or suite == "static1") {
    return {{autopas::OpenMPKindOption::omp_static}, {1}};
  }
  if (suite == "core") {
    return {{autopas::OpenMPKindOption::omp_static, autopas::OpenMPKindOption::omp_dynamic,
             autopas::OpenMPKindOption::lb4omp_static_steal},
            {1, autopas::OpenMPChunkSize::balanced, autopas::OpenMPChunkSize::golden}};
  }
  if (suite == "extended") {
    return {{autopas::OpenMPKindOption::omp_static, autopas::OpenMPKindOption::omp_dynamic,
             autopas::OpenMPKindOption::omp_guided, autopas::OpenMPKindOption::omp_auto,
             // LB4OMP's scheduling techniques used by Auto4OMP:
             autopas::OpenMPKindOption::llvm_trapezoidal, autopas::OpenMPKindOption::lb4omp_fac2a,
             autopas::OpenMPKindOption::lb4omp_static_steal, autopas::OpenMPKindOption::lb4omp_awf_b,
             autopas::OpenMPKindOption::lb4omp_awf_c, autopas::OpenMPKindOption::lb4omp_awf_d,
             autopas::OpenMPKindOption::lb4omp_awf_e, autopas::OpenMPKindOption::lb4omp_af_a},
            {1, 2, 1000, autopas::OpenMPChunkSize::balanced, autopas::OpenMPChunkSize::golden}};
  }
  autopas::utils::ExceptionHandler::exception(
      "getOpenMPSchedulesToTest(): Unknown value \"{}\" of AUTOPAS_TEST_OPENMP_SCHEDULES. Possible values: static1, "
      "core, extended.",
      suite);
  return {};
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
 * @param allowedOpenMPKindOptions By default, the kinds of getOpenMPSchedulesToTest(). Sequential traversals always use
 * the sequential kind, which is included implicitly.
 * @param allowedOpenMPChunkSizes By default, the chunk sizes of getOpenMPSchedulesToTest(). Sequential traversals
 * always use the not-applicable chunk size, which is included implicitly.
 * @param allowedVectorPatterns By default, all options.
 * @param throwIfNone If true (default), throw when no valid configuration exists; if false, return an empty set.
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
    const std::set<double> &allowedCellSizeFactors = {0.5, 1.0, 1.5},
    const std::set<autopas::OpenMPKindOption> &allowedOpenMPKindOptions = getOpenMPSchedulesToTest().kinds,
    const std::set<autopas::OpenMPChunkSize> &allowedOpenMPChunkSizes = getOpenMPSchedulesToTest().chunkSizes,
    const std::set<autopas::VectorizationPatternOption> &allowedVectorPatterns =
        autopas::VectorizationPatternOption::getAllOptions(),
    bool throwIfNone = true) {
  const autopas::NumberSetFinite<double> csfs(allowedCellSizeFactors);
  if (interactionType == autopas::InteractionTypeOption::all) {
    std::set<autopas::Configuration> allConfigs;
    for (auto iType : autopas::InteractionTypeOption::getMostOptions()) {
      // Do not throw for individual interaction types, as only the union has to be non-empty.
      const auto configs = autopas::SearchSpaceGenerators::cartesianProduct(
          allowedContainerOptions, allowedTraversalOptions, allowedLoadEstimatorOptions, allowedDataLayoutOptions,
          allowedNewton3Options, &csfs, allowedOpenMPKindOptions, allowedOpenMPChunkSizes, allowedVectorPatterns, iType,
          /*throwIfNone*/ false);
      allConfigs.insert(configs.begin(), configs.end());
    }
    if (throwIfNone and allConfigs.empty()) {
      autopas::utils::ExceptionHandler::exception("generateAllValidConfigurations: No valid configuration exists.");
    }
    return allConfigs;
  } else {
    return autopas::SearchSpaceGenerators::cartesianProduct(
        allowedContainerOptions, allowedTraversalOptions, allowedLoadEstimatorOptions, allowedDataLayoutOptions,
        allowedNewton3Options, &csfs, allowedOpenMPKindOptions, allowedOpenMPChunkSizes, allowedVectorPatterns,
        interactionType, throwIfNone);
  }
}

/**
 * Generates an arbitrary valid configuration constrained to the provided arguments. Where nullopt is given, no
 * constraints are applied.
 *
 * @param interactionType If provided, restrict to this interaction type, otherwise consider all interaction types.
 * @param containerOption If provided, restrict to this container, otherwise consider all options.
 * @param traversalOption If provided, restrict to this traversal, otherwise consider all options.
 * @param loadEstimatorOption If provided, restrict to this load estimator, otherwise consider all options.
 * @param dataLayoutOption If provided, restrict to this data layout, otherwise consider all options.
 * @param newton3Option If provided, restrict to this Newton 3 option, otherwise consider all options.
 * @param cellSizeFactor If provided, restrict to this cell size factor, otherwise consider {0.5, 1.0, 1.5}.
 * @param ompKind If provided, restrict to this OpenMP schedule kind, otherwise only static is considered. Sequential
 * traversals always use the sequential kind, regardless of this.
 * @param ompChunkSize If provided, restrict to this OpenMP chunk size, otherwise only 1 is considered. Sequential
 * traversals always use the not-applicable chunk size, regardless of this.
 * @param vectorPattern If provided, restrict to this vectorization pattern, otherwise consider all options.
 * @param throwIfNone If true (default), throw when no valid configuration exists; if false, return std::nullopt
 * instead.
 * @return An arbitrary valid configuration matching the given constraints, or std::nullopt if none exists.
 */
inline std::optional<autopas::Configuration> getArbitraryConfiguration(
    std::optional<autopas::InteractionTypeOption> interactionType = std::nullopt,
    std::optional<autopas::ContainerOption> containerOption = std::nullopt,
    std::optional<autopas::TraversalOption> traversalOption = std::nullopt,
    std::optional<autopas::LoadEstimatorOption> loadEstimatorOption = std::nullopt,
    std::optional<autopas::DataLayoutOption> dataLayoutOption = std::nullopt,
    std::optional<autopas::Newton3Option> newton3Option = std::nullopt,
    std::optional<double> cellSizeFactor = std::nullopt,
    std::optional<autopas::OpenMPKindOption> ompKind = std::nullopt,
    std::optional<autopas::OpenMPChunkSize> ompChunkSize = std::nullopt,
    std::optional<autopas::VectorizationPatternOption> vectorPattern = std::nullopt, bool throwIfNone = true) {
  const auto configurations = generateAllValidConfigurations(
      interactionType.value_or(autopas::InteractionTypeOption::all),
      containerOption.has_value() ? std::set<autopas::ContainerOption>{*containerOption}
                                  : autopas::ContainerOption::getAllOptions(),
      traversalOption.has_value() ? std::set<autopas::TraversalOption>{*traversalOption}
                                  : autopas::TraversalOption::getAllOptions(),
      loadEstimatorOption.has_value() ? std::set<autopas::LoadEstimatorOption>{*loadEstimatorOption}
                                      : autopas::LoadEstimatorOption::getAllOptions(),
      dataLayoutOption.has_value() ? std::set<autopas::DataLayoutOption>{*dataLayoutOption}
                                   : autopas::DataLayoutOption::getAllOptions(),
      newton3Option.has_value() ? std::set<autopas::Newton3Option>{*newton3Option}
                                : autopas::Newton3Option::getAllOptions(),
      cellSizeFactor.has_value() ? std::set<double>{*cellSizeFactor} : std::set<double>{0.5, 1.0, 1.5},
      ompKind.has_value() ? std::set<autopas::OpenMPKindOption>{*ompKind}
                          : std::set<autopas::OpenMPKindOption>{autopas::OpenMPKindOption::omp_static},
      ompChunkSize.has_value() ? std::set<autopas::OpenMPChunkSize>{*ompChunkSize}
                               : std::set<autopas::OpenMPChunkSize>{1},
      vectorPattern.has_value() ? std::set<autopas::VectorizationPatternOption>{*vectorPattern}
                                : autopas::VectorizationPatternOption::getAllOptions(),
      /*throwIfNone -> false so that we use the error message below which is clearer for users of this function*/
      false);

  if (configurations.empty()) {
    if (throwIfNone) {
      autopas::utils::ExceptionHandler::exception(
          "getArbitraryConfiguration: No valid configuration exists for the given constraints.");
    }
    return std::nullopt;
  }
  return *configurations.begin();
}

/**
 * Struct to hold a container and a cell size factor.
 */
struct ContainerConfiguration {
  /**
   * The container option.
   */
  autopas::ContainerOption container;
  /**
   * The cell size factor.
   */
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

  /**
   * Short string representation of the container configuration that is safe for test names.
   *
   * @note All "-" and "." are replaced with "_" so the result is safe to use as (part of) a test name.
   * @return
   */
  [[nodiscard]] std::string toShortString() const {
    std::string shortString = container.to_string() + "_csf_" + std::to_string(cellSizeFactor);
    std::ranges::replace(shortString, '-', '_');
    std::ranges::replace(shortString, '.', '_');
    return shortString;
  }

  /**
   * Generates an arbitrary valid full configuration using this container and cell size factor.
   *
   * @param interactionType If provided, restrict to this interaction type, otherwise consider all interaction types.
   * @param throwIfNone If true (default), throw when no valid configuration exists; if false, return std::nullopt.
   * @return An arbitrary valid full configuration using this container and cell size factor, or std::nullopt if none
   * exists..
   */
  [[nodiscard]] std::optional<autopas::Configuration> generateFullConfig(
      std::optional<autopas::InteractionTypeOption> interactionType = std::nullopt, bool throwIfNone = true) const {
    return getArbitraryConfiguration(interactionType, container, std::nullopt, std::nullopt, std::nullopt, std::nullopt,
                                     cellSizeFactor, /*ompKind*/ std::nullopt, /*ompChunkSize*/ std::nullopt,
                                     /*vectorPattern*/ std::nullopt, throwIfNone);
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
      const ContainerConfiguration containerConfig{containerOption, csf};
      // Keep this container configuration only if at least one valid full configuration exists for it.
      if (containerConfig.generateFullConfig(std::nullopt, /*throwIfNone*/ false).has_value()) {
        containerConfigs.insert(containerConfig);
      }
    }
  }
  return containerConfigs;
}
