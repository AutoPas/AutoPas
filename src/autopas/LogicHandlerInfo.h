/**
 * @file LogicHandlerInfo.h
 * @author F. Gratl
 * @date 19/06/2023
 */

#pragma once

#include "array"
#include "string"

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
  double verletSkin{0.4};
  /**
   * Number of particles in a cluster to use in VCL.
   */
  unsigned int verletClusterSize{4};
  /**
   * Number of particles in two cells from which sorting should be performed for traversal that use the CellFunctor
   */
  size_t sortingThreshold{8};
  /**
   * Time step used in the simulation.
   * This is currently used in rebuild frequency estimation for dynamic containers.
   */
  double deltaT{0};

  bool orderCellsByMortonIndex = false;
  bool useOptimizedLJFunctor = false;
  bool useCompactSoA = false;
  bool reserveVLSizes = false;
  bool bucketSortParticles = false;
  bool sortVerletLists = false;
  size_t sortingFrequency = 1;
};
}  // namespace autopas
