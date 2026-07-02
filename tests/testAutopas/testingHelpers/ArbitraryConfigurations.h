/**
 * @file ArbitraryConfigurations.h
 * @author S. Newcome
 * @date 02.07.2026
 */

#pragma once

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/tuning/Configuration.h"

/**
 * Arbitrary configurations for use in tests where the actual configuration does not matter, but where we need X unique
 * configurations.
 *
 * We distinguish still between pairwise and triwise configurations.
 *
 * @note The concrete option values of every configuration below are ARBITRARY. Tests using them should not assume that
 * the configuration will not change, other than whether it is pairwise or triwise.
 */
namespace arbitraryConfigurations {

// Pairwise configurations:
inline constexpr autopas::Configuration _arbitrary_config_2B_0 = autopas::Configuration(
    autopas::ContainerOption::directSum, 1.0, autopas::TraversalOption::ds_sequential,
    autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled,
    autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::NA);
inline constexpr autopas::Configuration _arbitrary_config_2B_1 = autopas::Configuration(
    autopas::ContainerOption::directSum, 1.0, autopas::TraversalOption::ds_sequential,
    autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
    autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::NA);
inline constexpr autopas::Configuration _arbitrary_config_2B_2 = autopas::Configuration(
    autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
    autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, autopas::InteractionTypeOption::pairwise,
    autopas::VectorizationPatternOption::NA);
inline constexpr autopas::Configuration _arbitrary_config_2B_3 = autopas::Configuration(
    autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
    autopas::VectorizationPatternOption::NA);
inline constexpr autopas::Configuration _arbitrary_config_2B_4 = autopas::Configuration(
    autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c18, autopas::LoadEstimatorOption::none,
    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
    autopas::VectorizationPatternOption::NA);
inline constexpr autopas::Configuration _arbitrary_config_2B_5 = autopas::Configuration(
    autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c18, autopas::LoadEstimatorOption::none,
    autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, autopas::InteractionTypeOption::pairwise,
    autopas::VectorizationPatternOption::NA);
inline constexpr autopas::Configuration _arbitrary_config_2B_6 = autopas::Configuration(
    autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
    autopas::VectorizationPatternOption::NA);

// Triwise (three-body) configurations:
inline constexpr autopas::Configuration _arbitrary_config_3B_0 = autopas::Configuration(
    autopas::ContainerOption::linkedCells, 1.0, autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::triwise,
    autopas::VectorizationPatternOption::NA);
inline constexpr autopas::Configuration _arbitrary_config_3B_1 = autopas::Configuration(
    autopas::ContainerOption::directSum, 1.0, autopas::TraversalOption::ds_sequential,
    autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
    autopas::InteractionTypeOption::triwise, autopas::VectorizationPatternOption::NA);

}  // namespace arbitraryConfigurations