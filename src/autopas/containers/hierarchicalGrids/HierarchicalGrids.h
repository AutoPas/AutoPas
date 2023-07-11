/**
 * @file HierarchicalGrids.h
 * @author hoppef
 * @brief 
 * @version 0.1
 * @date 2023-07-10
 * 
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCells.h"

#include <vector>

namespace autopas {

template <class Particle>
class HierarchicalGrids : public ParticleContainerInterface<FullParticleCell<Particle>> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = typename ParticleCell::ParticleType;

  HierarchicalGrids (const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
                      const double skinPerTimestep, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0,
                      LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell, const unsigned int numberOfLevels)
      : _numberOfLevels(numberOfLevels) {

    // Loop over all hierarchy levels and create a Linked Cells container for each level.
    //
    // @todo !!! cellSizeFactor has to be adapted !!!
    //
    hierarchyLevels.reserve(_numberOfLevels);
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels.emplace_back(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, cellSizeFactor,loadEstimator);
    }
    
  }

  /**
   * @todo !
   * @brief Get the ContainerType.
   * @return ContainerOption of the type of this container.
   */
  [[nodiscard]] virtual ContainerOption getContainerType() const = 0;

    /**
   *
   * Get the ParticleCell type as an Enum
   * @return The Cell type as an Enum
   */
  [[nodiscard]] CellType getParticleCellTypeEnum() override { return CellType::FullParticleCell; }

    /**
   * @todo !
   * @copydoc AutoPas::reserve()
   * @param numParticlesHaloEstimate Estimate for the number of halo particles.
   * Reserves space in the container data structure.
   */
  virtual void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {

    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels[level].reserve(numParticles, numParticlesHaloEstimate);
    }

  }

  /**
   * 
   * @todo !
   * 
   * Adds a particle to the container.
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param p The particle to be added. This particle is already checked to be inside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be inside of the bounding box!
   */
  void addParticleImpl(const Particle &p) override {

    // @todo !!! Write getter !!!
    unsigned int level = getHierarchyLevelOfParticle(p);
    hierarchyLevels[level].addParticle(p);

  }

  /**
   * 
   * @todo !
   * 
   * Adds a particle to the container that lies in the halo region of the container.
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param haloParticle Particle to be added. This particle is already checked to be outside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be outside of the bounding box!
   */
  void addHaloParticleImpl(const Particle &haloParticle) override {

    // @todo !!! Write getter !!!
    unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    hierarchyLevels[level].addHaloParticle(haloParticle);

  }

    /**
   *    
   * @todo !
   * 
   * Update a halo particle of the container with the given haloParticle.
   * @param haloParticle Particle to be updated.
   * @return Returns true if the particle was updated, false if no particle could be found.
   */
  virtual bool updateHaloParticle(const Particle &haloParticle) = 0;

  /**
   * @todo !
   * Rebuilds the neighbor lists.
   * @param traversal The used traversal.
   */
  virtual void rebuildNeighborLists(TraversalInterface *traversal) = 0;

  /**
   * 
   * @todo !
   * Deletes all halo particles.
   */
  virtual void deleteHaloParticles() = 0;

  /**
   * @todo !
   * Deletes all particles.
   */
  virtual void deleteAllParticles() = 0;

  /**
   * @todo !
   * Get the total number of particles saved in the container (owned and halo but not dummies).
   * @return Number of particles in the container.
   */
  [[nodiscard]] virtual unsigned long getNumberOfParticles() const = 0;



  private:

/**
 * @brief Vector containing the HGrid's different hierarchy levels of Linked Cells.
 * 
 */
  std::vector<LinkedCells<ParticleCell>> hierarchyLevels;

  /**
   * @brief Number of levels (hierarchies) in the H-Grid
   * 
   */
  unsigned int _numberOfLevels = 1;

};
  
} // namespace autopas 