/**
 * @file SearchSpaceGenerators.h
 * @author F. Gratl
 * @date 23.06.23
 */

#pragma once

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/utils//Evidence.h"

namespace autopas::SearchSpaceGenerators {
/**
 * Fills the search space with the cartesian product of the given options (minus invalid combinations).
 * @param allowedContainerOptions
 * @param allowedCellSizeFactors
 * @param allowedTraversalOptions
 * @param allowedLoadEstimatorOptions
 * @param allowedDataLayoutOptions
 * @param allowedNewton3Options
 * @return A map containing all valid configurations pointing to empty vectors. Map<Configuration, {}>
 */
std::map<autopas::Configuration, std::vector<autopas::Evidence>> optionCrossProduct(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
    const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options);

}  // namespace autopas::SearchSpaceGenerators
