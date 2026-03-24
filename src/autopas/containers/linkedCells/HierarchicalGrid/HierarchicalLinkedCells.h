/**
 * @file HierarchicalLinkedCells.h
 *
 * @date 21 Mar 2026
 */

#pragma once

#include "autopas/containers/linkedCells/HierarchicalGrid/HierarchicalCellBlock3D.h"
#include "autopas/containers/linkedCells/LinkedCells.h"

namespace autopas {

/**
 * LinkedCells variant that uses HierarchicalCellBlock3D.
 * @tparam Particle_T type of the Particle
 */
template <class Particle_T>
class HierarchicalLinkedCells : public LinkedCells<Particle_T> {
 public:
  /**
   * Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle_T>;

  /**
   * Type of the Particle.
   */
  using ParticleType = typename ParticleCell::ParticleType;

  /**
   * Constructor of the HierarchicalLinkedCells class.
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param rebuildFrequency
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param loadEstimator the load estimation algorithm for balanced traversals.
   */
  HierarchicalLinkedCells(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
              const double skin, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0,
              LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell)
      : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skin, rebuildFrequency),
        this->_cellBlock(this->_cells, boxMin, boxMax, cutoff + skin, cellSizeFactor),
        this->_loadEstimator(loadEstimator) {}
};

}  // namespace autopas
