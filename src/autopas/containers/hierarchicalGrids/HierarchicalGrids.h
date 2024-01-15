/**
 * @file HierarchicalGrids.h
 * @author hoppef
 * @date 2023-07-10
 * 
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/options/LoadEstimatorOption.h"

#include "../applicationLibrary/discreteElementMethod/discreteElementMethodLibrary/DEMFunctor.h"

#include <vector>
#include <deque>

namespace autopas {

template <class Particle>
class HierarchicalGrids : public ParticleContainerInterface<Particle> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = Particle; //FullParticleCell<Particle>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle; // debugging mit Piet
  //using ParticleType = typename ParticleCell::ParticleType;

  HierarchicalGrids (const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
                      const double skinPerTimestep, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0,
                      LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell, const unsigned int numberOfLevels = 1)
      : _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff), _skinPerTimestep(skinPerTimestep), _rebuildFrequency(rebuildFrequency), _cellSizeFactor(cellSizeFactor), _loadEstimator(loadEstimator), _numberOfLevels(numberOfLevels) {

    // Loop over all hierarchy levels and create a Linked Cells container for each level.
    //
    // @todo !!! cellSizeFactor has to be adapted !!!
    //
    //hierarchyLevels.reserve(_numberOfLevels);
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels.emplace_back(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, cellSizeFactor, loadEstimator);
    }
    
  }

  /**
   * @brief Get the ContainerType.
   * @return ContainerOption of the type of this container.
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::hierarchicalGrids; };

    /**
   * Get the ParticleCell type as an Enum
   * @return The Cell type as an Enum
   */
  [[nodiscard]] CellType getParticleCellTypeEnum() const override { return CellType::FullParticleCell; }

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

  protected:
  /**
   * @brief Adds a particle to the container.
   * 
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param p The particle to be added. This particle is already checked to be inside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be inside of the bounding box!
   */
  void addParticleImpl(const ParticleType &p) override {

    //unsigned int level = getHierarchyLevelOfParticle(p);
    //addParticleToGridLevel(p, level);

    hierarchyLevels[0].addParticleImpl(p);

  }

  public:
  /**
   * @brief Adds a particle to a given level of the hierarchical grid
   * 
   * @param p The particle to be added.
   * @param level The level of the H-grid to be added to.
   */
  void addParticleToGridLevel(const ParticleType &p, unsigned int level) {

    hierarchyLevels[level].addParticle(p);

  }

  protected:
  /**
   * @brief Adds a particle to the container that lies in the halo region of the container.
   * 
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param haloParticle Particle to be added. This particle is already checked to be outside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be outside of the bounding box!
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {

    //unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    //hierarchyLevels[level].addHaloParticle(haloParticle);

    hierarchyLevels[0].addHaloParticle(haloParticle);

  }

  public:
  /** 
   * @brief Update a halo particle of the container with the given haloParticle.
   * 
   * @param haloParticle Particle to be updated.
   * @return Returns true if the particle was updated, false if no particle could be found.
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {

    //unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    //return hierarchyLevels[level].updateHaloParticle(haloParticle);

    return hierarchyLevels[0].updateHaloParticle(haloParticle);

  }

  /**
   * Rebuilds the neighbor lists.
   * Nothing to do here.
   * @param traversal The used traversal.
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
  }

  /**
   * Deletes all halo particles.
   */
  void deleteHaloParticles() override {
    //for (unsigned int level = 0; level < _numberOfLevels; level++) {
    //  hierarchyLevels[level].deleteHaloParticles();
    //}
    hierarchyLevels[0].deleteHaloParticles();
  }

  /**
   * Deletes all particles.
   */
  void deleteAllParticles() override {
    //for (unsigned int level = 0; level < _numberOfLevels; level++) {
    //  hierarchyLevels[level].deleteAllParticles();
    //}
    hierarchyLevels[0].deleteAllParticles();
  }

  /**
   * Get the number of particles with respect to the specified IteratorBehavior.
   * @warning: Since this function counts the number of the respective particles in the internal particle storage, this
   * is in O(n) + lock is required. Only use it when it is absolutely necessary to have the exact number of different
   * particle types like owned or halo. If it is enough to have the whole number of particles (owned + halo + dummy),
   * the function size() can be used.
   * @param behavior Behavior of the iterator, see IteratorBehavior.
   * @return Number of particles in the container.
   */
  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const override {
    //size_t numParticles = 0ul;
    //for (unsigned int level = 0; level < _numberOfLevels; level++) {
    //  numParticles += hierarchyLevels[level].getNumberOfParticles(behavior);
    //}
    //return numParticles;
    return hierarchyLevels[0].getNumberOfParticles(behavior);
    
  }

  /**
   * Get the total number of particles saved in the container (owned + halo + dummy).
   * @return Number of particles saved in the container (owned + halo + dummy).
   */
  [[nodiscard]] size_t size() const override {
    size_t numParticles = 0ul;
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      numParticles += hierarchyLevels[level].size();
    }
    return numParticles;
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<ParticleType, true, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, true, false>::ParticleVecType *additionalVectors = nullptr) override {
        
        return hierarchyLevels[0].begin(behavior, additionalVectors);

      }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note const version
   */
  [[nodiscard]] ContainerIterator<ParticleType, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, false, false>::ParticleVecType *additionalVectors = nullptr) const override {
        
        return hierarchyLevels[0].begin(behavior, additionalVectors);

      }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<ParticleType, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, true, true>::ParticleVecType *additionalVectors = nullptr) override {

        return hierarchyLevels[0].getRegionIterator(lowerCorner, higherCorner, behavior, additionalVectors);

      }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   * @note const version
   */
  [[nodiscard]] ContainerIterator<ParticleType, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, false, true>::ParticleVecType *additionalVectors = nullptr) const override {

        return hierarchyLevels[0].getRegionIterator(lowerCorner, higherCorner, behavior, additionalVectors);

      }

  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    
    hierarchyLevels[0].forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);

    }

  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    
    hierarchyLevels[0].reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);

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
    unsigned int numberOfCrossLevels = getNumberOfCrossLevels(_numberOfLevels);

    //@todo cellSizeFactor has to be adapted

    for (unsigned int largerLevel = 0; largerLevel < _numberOfLevels - 1; largerLevel++) {
      for (unsigned int smallerLevel = largerLevel; smallerLevel < _numberOfLevels; smallerLevel++) {
        
        //Create cross-level
        autopas::LinkedCells<ParticleCell> crossLevel(_boxMin, _boxMax, _cutoff, _skinPerTimestep, _rebuildFrequency, _cellSizeFactor, _loadEstimator);
        
        //Add large particles
        for(auto iterLarger = hierarchyLevels[largerLevel].begin(); iterLarger.isValid(); ++iterLarger) {
          crossLevel.addParticle(*iterLarger);
        }

        //Add small particles
        for(auto iterSmaller = hierarchyLevels[smallerLevel].begin(); iterSmaller.isValid(); ++iterSmaller) {
          crossLevel.addParticle(*iterSmaller);
        }

        // Iterate pairwise over the cross-level
        crossLevel.iteratePairwise(traversal);

      }
    }
    

    // Iterate over hierarchy levels and subtract excess forces
    demLib::DEMFunctor<Particle>::initExcessForceSubtraction(_numberOfLevels);
    for (unsigned int level = 0; level < _numberOfLevels; level++) {
      hierarchyLevels[level].iteratePairwise(traversal); // add traversal option
    }
    demLib::DEMFunctor<Particle>::endExcessForceSubtraction();

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMax()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const final { return _boxMax; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMin()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const final { return _boxMin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getCutoff()
   */
  [[nodiscard]] double getCutoff() const final { return _cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) final { _cutoff = cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getVerletSkin()
   */
  [[nodiscard]] double getVerletSkin() const override { return _skinPerTimestep * _rebuildFrequency; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
   */
  [[nodiscard]] double getInteractionLength() const override { return _cutoff + _skinPerTimestep * _rebuildFrequency; }

  /**
   * @TODO !!! Funktioniert so bestimmt nicht ... 
   * 
   * @copydoc autopas::ParticleContainerInterface::updateContainer()
   */
  [[nodiscard]] std::vector<ParticleType> updateContainer(bool keepNeighborListsValid) {

    return hierarchyLevels[0].updateContainer(keepNeighborListsValid);

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {

    return hierarchyLevels[0].getTraversalSelectorInfo();

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getParticle()
   */
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                            IteratorBehavior iteratorBehavior) const override {

    return hierarchyLevels[0].getParticle(cellIndex, particleIndex, iteratorBehavior); 

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getParticle()
   */
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                                   IteratorBehavior iteratorBehavior,
                                                                   const std::array<double, 3> &boxMin,
                                                                   const std::array<double, 3> &boxMax) const override {

    return hierarchyLevels[0].getParticle(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::deleteParticle()
   */
  bool deleteParticle(Particle &particle) override {

    return hierarchyLevels[0].deleteParticle(particle);

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::deleteParticle()
   */
  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {

    return hierarchyLevels[0].deleteParticle(cellIndex, particleIndex);

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

    unsigned int level;

    // Loop over all hierarchy levels
    for (level = 0; level < _numberOfLevels; level++) {
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
  //static std::vector<LinkedCells<Particle>> hierarchyLevels;
  std::deque<LinkedCells<Particle>> hierarchyLevels;

  //std::vector<LinkedCells<Particle>> crossLevels;

  /**
   * @brief Number of levels (hierarchies) in the H-Grid
   * 
   */
  unsigned int _numberOfLevels = 1;

  std::array<double, 3> _boxMin; 
  std::array<double, 3> _boxMax; 
  double _cutoff;
  double _skinPerTimestep; 
  unsigned int _rebuildFrequency;
  double _cellSizeFactor;
  autopas::LoadEstimatorOption _loadEstimator;

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