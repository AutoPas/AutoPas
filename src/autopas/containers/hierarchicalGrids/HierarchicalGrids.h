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
      hierarchyLevels.emplace_back(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, cellSizeFactor, loadEstimator);
    }
    
  }

  /**
   * 
   * @brief Get the ContainerType.
   * @return ContainerOption of the type of this container.
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::hierarchicalGrids; };

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
  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {

    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels[level].reserve(numParticles, numParticlesHaloEstimate);
    }

  }

  /**
   * @brief Adds a particle to the container.
   * 
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param p The particle to be added. This particle is already checked to be inside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be inside of the bounding box!
   */
  void addParticleImpl(const ParticleType &p) override {

    unsigned int level = getHierarchyLevelOfParticle(p);
    addParticleToGridLevel(p, level);

  }

  /**
   * @brief Adds a particle to a given level of the hierarchical grid
   * 
   * @param p The particle to be added.
   * @param level The level of the H-grid to be added to.
   */
  void addParticleToGridLevel(const ParticleType &p, unsigned int level) {

    hierarchyLevels[level].addParticle(p);

  }

  /**
   * 
   * @brief Adds a particle to the container that lies in the halo region of the container.
   * 
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param haloParticle Particle to be added. This particle is already checked to be outside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be outside of the bounding box!
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {

    unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    hierarchyLevels[level].addHaloParticle(haloParticle);

  }

  /**
   *    
   * @brief Update a halo particle of the container with the given haloParticle.
   * 
   * @param haloParticle Particle to be updated.
   * @return Returns true if the particle was updated, false if no particle could be found.
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {

    unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    hierarchyLevels[level].updateHaloParticle(haloParticle);

  }

  /**
   * 
   * Rebuilds the neighbor lists.
   * Nothing to do here.
   * @param traversal The used traversal.
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
  }
  /**
   * 
   * Deletes all halo particles.
   */
  void deleteHaloParticles() override {
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels[level].deleteHaloParticles();
    }
  }

  /**
   * 
   * Deletes all particles.
   */
  void deleteAllParticles() override {
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels[level].deleteAllParticles();
    }
  }

  /**
   * @todo !
   * Get the total number of particles saved in the container (owned and halo but not dummies).
   * @return Number of particles in the container.
   */
  [[nodiscard]] unsigned long getNumberOfParticles() override {
    size_t numParticles = 0ul;
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      numParticles += hierarchyLevels[level].getNumberOfParticles();
    }
    return numParticles;
  }

  /**
   * @todo
   * 
   * Iterates over all particle pairs in the container.
   * @param traversal The traversal to use for the iteration.
   */
  void iteratePairwise(TraversalInterface *traversal) override {

    // Iterate over cross levels

    // Get the number of hierarchy cross-levels
    numberOfCrossLevels = getNumberOfCrossLevels(_numberOfLevels);

    //@todo cellSizeFactor has to be adapted
    /*
    crossLevels.reserve(numberOfCrossLevels);
    for (unsigned int level = 0; level < numberOfCrossLevels; level++) {
      crossLevels.emplace_back(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, cellSizeFactor, loadEstimator);
    }
    */

    for (unsigned int largerLevel = 0; largerLevel < _numberOfLevels - 1; largerLevel++) {
      for (unsigned int smallerLevel = largerLevel; smallerLevel < _numberOfLevels; smallerLevel++) {
        
        //Create cross-level
        autopas::LinkedCells<ParticleCell> crossLevel(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, cellSizeFactor, loadEstimator);
        
        //Add large particles
        for(auto iterLarger = hierarchyLevels[largerLevel].begin(); iterLarger.isValid(); ++iterLarger) {
          crossLevel.addParticle(*iterLarger);
        }

        //Add small particles
        for(auto iterSmaller = hierarchyLevels[smallerLevel].begin(); iterSmaller.isValid(); ++iterSmaller) {
          crossLevel.addParticle(*iterSmaller);
        }

        crossLevel.iteratePairwise(&traversal);

      }
    }
    

    // Iterate over hierarchy levels and subtract excess forces
    autopas::DEMFunctor::initExcessForceSubtraction(_numberOfLevels);
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels[level].iteratePairwise(&traversal); // add traversal option
    }
    autopas::DEMFunctor::endExcessForceSubtraction();

  }

  /**
   * @brief Get the Hierarchy Level Of Particle object
   * 
   * @param p 
   * @return Level the particle should be sorted into. 
   */
  unsigned int getHierarchyLevelOfParticle(const ParticleType &p) {
    
    // Get the radius of the particle in question
    auto particleRad = p.getRad();

    // Loop over all hierarchy levels
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      // As long as the 'hierarchicalGridBorders' element is smaller than the radius continue to cycle
      if (particleRad > hierarchicalGridBorders[level]) break;     
    }

    if(level == 0) {
    utils::ExceptionHandler::exception(
      "HierarchicalGrids: Trying to get the hierarchy level of a particle which is too large to fit in the grid.\n"
      "Particle radius {}\n"
      "Upper limit of the largest level {}",
      particleRad, hierarchicalGridBorders[level]);
    }
    // Break and return the level
    return level - 1;
  }


  private:

  /**
   * @brief Get the number Of cross-levels of the H-Grid
   * 
   * @param levels 
   * @return unsigned int 
   */
  unsigned int getNumberOfCrossLevels(unsigned int levels) {

    unsigned int numberOfCrossLevels = 0;
    for (unsigned int i = 1; i < levels; i++)
    {
      numberOfCrossLevels += i;
    }
    
    return numberOfCrossLevels;
  }

  /**
  * @brief Vector containing the HGrid's different hierarchy levels of Linked Cells.
  * 
  */
  std::vector<LinkedCells<ParticleCell>> hierarchyLevels;

  std::vector<LinkedCells<ParticleCell>> crossLevels;

  /**
   * @brief Number of levels (hierarchies) in the H-Grid
   * 
   */
  unsigned int _numberOfLevels = 1;

  /**
   * @brief Vector containing the upper boundaries of the grid levels. 
   * 
   * For example, the first element describes the upper boundary of the largest H-Grid level,
   * i.e. the largest particles fitting in the topmost level.
   * By design 'hierarchicalGridBorders' should have the same '_numberOfLevels' as the H-Grid.
   * 
   */
  std::vector<double> hierarchicalGridBorders;

};
  
} // namespace autopas 