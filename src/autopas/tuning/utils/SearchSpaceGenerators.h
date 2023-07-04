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
#include "autopas/utils/NumberSetFinite.h"

namespace autopas::SearchSpaceGenerators {

/**
 * Helper struct for all dimensions of a full cross product search space.
 */
struct OptionSpace {
  /*
   * Available container options.
   */
  std::set<ContainerOption> containerOptions;
  /*
   * Available traversal options.
   */
  std::set<TraversalOption> traversalOptions;
  /*
   * Available loadEstimator options.
   */
  std::set<LoadEstimatorOption> loadEstimatorOptions;
  /*
   * Available dataLayout options.
   */
  std::set<DataLayoutOption> dataLayoutOptions;
  /*
   * Available newton3 options.
   */
  std::set<Newton3Option> newton3Options;
  /*
   * Available discrete cellSizeFactors options.
   */
  std::set<double> cellSizeFactors;
};

/**
 * Fills the search space with the cartesian product of the given options (minus invalid combinations).
 * @param allowedContainerOptions
 * @param allowedTraversalOptions
 * @param allowedLoadEstimatorOptions
 * @param allowedDataLayoutOptions
 * @param allowedNewton3Options
 * @param allowedCellSizeFactors
 * @return A set containing all valid configurations.
 */
std::set<Configuration> optionCrossProduct(const std::set<ContainerOption> &allowedContainerOptions,
                                           const std::set<TraversalOption> &allowedTraversalOptions,
                                           const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                           const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                           const std::set<Newton3Option> &allowedNewton3Options,
                                           std::unique_ptr<NumberSet<double>> allowedCellSizeFactors);

/**
 * Crudely trying to reconstruct the dimensions of the search space from a given set of options.
 *
 * @note It is assumed that searchSet is a full cross product.
 *
 * @param searchSet
 * @return
 */
OptionSpace inferOptionDimensions(const std::set<Configuration> &searchSet);

/**
 * For a given domain parametrization, calculate which cell size factors (csf) in an interval actually are useful to
 * generate distinct numbers of cells.
 *
 * A cell is assumed to be of side length interactionLength * csf.
 *
 * @param numberInterval
 * @param interactionLength
 * @param domainLengthX Only consider the domain length in one dimension, since the csf is also one dimensional.
 * @return Set of cell size factors that yield different numbers of cells.
 */
std::set<double> calculateRelevantCsfs(const NumberInterval<double> &numberInterval, double interactionLength,
                                       double domainLengthX);

}  // namespace autopas::SearchSpaceGenerators
