/**
 * @file LogicHandlerInfo.h
 * @author F. Gratl
 * @date 19/06/2023
 */

#pragma once

#include "array"
#include "string"
#include "vector"

namespace autopas {
/**
 * Class that wraps all arguments for the logic handler to provide a more stable API.
 */
class LogicHandlerInfo {
 public:
  /**
   * Lower corner of the container without halo.
   */
  std::array<double, 3> boxMin{0., 0., 0.};
  /**
   * Upper corner of the container without halo.
   */
  std::array<double, 3> boxMax{0., 0., 0.};
  /**
   * Cutoff radius to be used in this simulation.
   */
  double cutoff{1.};
  /**
   * Length added to the cutoff for the Verlet lists' skin.
   */
  double verletSkinPerTimestep{0.02};
  /**
   * Number of particles in a cluster to use in VCL.
   */
  unsigned int verletClusterSize{4};
  /**
   * Number of particles in two cells from which sorting should be performed for traversal that use the CellFunctor
   */
  size_t sortingThreshold{8};
  /**
   * Cutoff radius for each level of a HierarchicalGrid
   */
  std::vector<double> cutoffs{};
};
}  // namespace autopas
