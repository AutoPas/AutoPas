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
                      LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell, 
                      const unsigned int numberOfLevels = 1)
      : _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff), _skinPerTimestep(skinPerTimestep), _rebuildFrequency(rebuildFrequency), 
        _cellSizeFactor(cellSizeFactor), _loadEstimator(loadEstimator) {

    // Loop over all hierarchy levels and create a Linked Cells container for each level.
    //
    // @todo !!! cellSizeFactor has to be adapted !!!
    //


      // Simplistic first shot at the size of the individual H-Grid levels
      //
      // @todo !!! Implment a more dynamic approach accounting for the number of particles in the system and 
      // allocate the levels according to their size and distribution
      //
      const double radiusMax = cutoff / 2.;
      const double radiusSegment = radiusMax / numberOfLevels;

      std::cout << "Created Hierarchical Grid with " << numberOfLevels << " levels" << std::endl;

    // Reserve appropriate number of elements for _hierarchicalGridBorders
    // Note: _hierarchyLevels has to be a std::deque, which does not support reserve()
    _hierarchicalGridBorders.reserve(numberOfLevels);

    for (unsigned int level = 0; level < numberOfLevels; level++) {

      const double gridBorder = radiusMax - radiusSegment * level;
      const double cellSize = gridBorder / radiusMax;
      _hierarchicalGridBorders.emplace_back(gridBorder);

      _hierarchyLevels.emplace_back(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, 
                                   cellSizeFactor, loadEstimator);

      std::cout << "Created Hierarchical Grid level: " << level << std::endl;
      std::cout << "  with max. allowed particle radii of " << _hierarchicalGridBorders[level] << std::endl;
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

    for (unsigned int level = 0; level < _hierarchyLevels.size(); level++) {
      _hierarchyLevels[level].reserve(numParticles, numParticlesHaloEstimate);
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

    _hierarchyLevels[0].addParticleImpl(p);

  }

  public:
  /**
   * @brief Adds a particle to a given level of the hierarchical grid
   * 
   * This relies on an unsafe version of addParticle() and does not perform a boundary check.
   * @param p The particle to be added.
   * @param level The level of the H-grid to be added to.
   */
  void addParticleToGridLevel(const ParticleType &p, unsigned int level) {

    _hierarchyLevels[level].addParticleImpl(p);

  }

    /**
   * @brief Adds a particle to a given level of the hierarchical grid
   * 
   * This relies on an unsafe version of addHaloParticle() and does not perform a boundary check.
   * @param p The particle to be added.
   * @param level The level of the H-grid to be added to.
   */
  void addHaloParticleToGridLevel(const ParticleType &p, unsigned int level) {

    _hierarchyLevels[level].addHaloParticleImpl(p);

  }

  protected:
  /**
   * @brief Adds a particle to the container that lies in the halo region of the container.
   * 
   * This is an unsafe version of addHaloParticle() and does not perform a boundary check.
   * @param haloParticle Particle to be added. This particle is already checked to be outside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be outside of the bounding box!
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {

    //unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    //_hierarchyLevels[level].addHaloParticle(haloParticle);

    _hierarchyLevels[0].addHaloParticleImpl(haloParticle);

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
    //return _hierarchyLevels[level].updateHaloParticle(haloParticle);

    return _hierarchyLevels[0].updateHaloParticle(haloParticle);

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
    //for (unsigned int level = 0; level < _hierarchyLevels.size(); level++) {
    //  _hierarchyLevels[level].deleteHaloParticles();
    //}
    _hierarchyLevels[0].deleteHaloParticles();
  }

  /**
   * Deletes all particles.
   */
  void deleteAllParticles() override {
    //for (unsigned int level = 0; level < _hierarchyLevels.size(); level++) {
    //  _hierarchyLevels[level].deleteAllParticles();
    //}
    _hierarchyLevels[0].deleteAllParticles();
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
    //for (unsigned int level = 0; level < _hierarchyLevels.size(); level++) {
    //  numParticles += _hierarchyLevels[level].getNumberOfParticles(behavior);
    //}
    //return numParticles;
    return _hierarchyLevels[0].getNumberOfParticles(behavior);
    
  }

  /**
   * Get the total number of particles saved in the container (owned + halo + dummy).
   * @return Number of particles saved in the container (owned + halo + dummy).
   */
  [[nodiscard]] size_t size() const override {
    /*size_t numParticles = 0ul;
    for (unsigned int level = 0; level < _hierarchyLevels.size(); level++) {
      numParticles += _hierarchyLevels[level].size();
    }
    return numParticles;*/
    return _hierarchyLevels[0].size();
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<ParticleType, true, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, true, false>::ParticleVecType *additionalVectors = nullptr) override {
        
        return _hierarchyLevels[0].begin(behavior, additionalVectors);

      }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note const version
   */
  [[nodiscard]] ContainerIterator<ParticleType, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, false, false>::ParticleVecType *additionalVectors = nullptr) const override {
        
        return _hierarchyLevels[0].begin(behavior, additionalVectors);

      }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<ParticleType, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, true, true>::ParticleVecType *additionalVectors = nullptr) override {

        return _hierarchyLevels[0].getRegionIterator(lowerCorner, higherCorner, behavior, additionalVectors);

      }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   * @note const version
   */
  [[nodiscard]] ContainerIterator<ParticleType, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, false, true>::ParticleVecType *additionalVectors = nullptr) const override {

        return _hierarchyLevels[0].getRegionIterator(lowerCorner, higherCorner, behavior, additionalVectors);

      }

  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    
    _hierarchyLevels[0].forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);

    }

  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    
    _hierarchyLevels[0].reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);

    }

    /**
   * Container specific implementation for getParticle. See ParticleContainerInterface::getParticle().
   *
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin
   * @param boxMax
   * @return tuple<ParticlePointer, CellIndex, ParticleIndex>
   */
  template <bool regionIter>
  std::tuple<const Particle *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
                                                               IteratorBehavior iteratorBehavior,
                                                               const std::array<double, 3> &boxMin,
                                                               const std::array<double, 3> &boxMax) const {
    // first and last relevant cell index
    const auto [startCellIndex, endCellIndex] = [&]() -> std::tuple<size_t, size_t> {
      if constexpr (regionIter) {
        return {_cellBlock.get1DIndexOfPosition(boxMin), _cellBlock.get1DIndexOfPosition(boxMax)};
      } else {
        if (not(iteratorBehavior & IteratorBehavior::halo)) {
          // only potentially owned region
          return {_cellBlock.getFirstOwnedCellIndex(), _cellBlock.getLastOwnedCellIndex()};
        } else {
          // whole range of cells
          return {0, this->_cells.size() - 1};
        }
      }
    }();

    // if we are at the start of an iteration ...
    if (cellIndex == 0 and particleIndex == 0) {
      cellIndex =
          startCellIndex + ((iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num());
    }
    // abort if the start index is already out of bounds
    if (cellIndex >= this->_cells.size()) {
      return {nullptr, 0, 0};
    }
    // check the data behind the indices
    if (particleIndex >= this->_cells[cellIndex].size() or
        not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
            this->_cells[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax)) {
      // either advance them to something interesting or invalidate them.
      std::tie(cellIndex, particleIndex) =
          advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax, endCellIndex);
    }

    // shortcut if the given index doesn't exist
    if (cellIndex > endCellIndex) {
      return {nullptr, 0, 0};
    }
    const Particle *retPtr = &this->_cells[cellIndex][particleIndex];

    return {retPtr, cellIndex, particleIndex};
  }



  void distributeParticlesToHGLevels() {

      AUTOPAS_OPENMP(parallel)
      for (auto iterParticle = _hierarchyLevels[0].begin(); iterParticle.isValid(); ++iterParticle) {
        // get level the particle should be inserted to
        unsigned int targetLevel = getHierarchyLevelOfParticle(*iterParticle);

        // No need to distribute particles, which are supposed to be in level = 0
        if (targetLevel == 0) { continue; }

        // add particles to their respective levels
        addParticleToGridLevel(*iterParticle, targetLevel);

        // and delete it from the carrier level
        internal::deleteParticle(iterParticle);
      }

  }

  void distributeParticlesToCarrierLevel() {

    //iterate over all HG levels
    for (unsigned int iterLevels = 1; iterLevels < _hierarchyLevels.size(); iterLevels++) {
      
      //iterate over all particles in a given level
      for(auto iterParticle = _hierarchyLevels[iterLevels].begin(); iterParticle.isValid(); ++iterParticle) {

        //add the particle to the carrier level
        _hierarchyLevels[0].addParticleImpl(*iterParticle);

      }
      _hierarchyLevels[iterLevels].deleteAllParticles();

    }

  }

  void iterateAllHGLevels(TraversalInterface *traversal) {

    // @idea: Reduce the number of threads per level and parallelize as much as 
    // possible over the levels. Might mean some major redesign all over autopas? 
    //AUTOPAS_OPENMP(parallel firstprivate (traversal))
    for (auto & hgLevel : _hierarchyLevels) {

      hgLevel.iteratePairwise(traversal);
      
    }
  }

  void iterateAllCrossLevels(TraversalInterface *traversal) {

    demLib::DEMFunctor<Particle>::initCrossLevelTraversal();

    for (unsigned int largerLevel = 0; largerLevel < _hierarchyLevels.size(); largerLevel++) {
      for (unsigned int smallerLevel = largerLevel + 1; smallerLevel < _hierarchyLevels.size(); smallerLevel++) {

        // If one of the HG levels is empty there is no need to traverse over the cross-level, 
        // since no particle interaction should be found. Hence, continue;
        const bool largerLevelEmpty = _hierarchyLevels[largerLevel].size() == 0;
        const bool smallerLevelEmpty = _hierarchyLevels[smallerLevel].size() == 0;
        if (largerLevelEmpty or smallerLevelEmpty){
          continue;
        }
        
        
        const double crossLevelCellSize = ( _hierarchyLevels[largerLevel].getCutoff() + 
                                            _hierarchyLevels[smallerLevel].getCutoff() ) / 2.0;

        //Create cross-level
        autopas::LinkedCells<Particle> crossLevel(_boxMin, _boxMax, crossLevelCellSize, _skinPerTimestep, _rebuildFrequency, 
                                                  _cellSizeFactor, _loadEstimator);
        
        //Very ugly trick: Make particle radii of the larger level negative
        //Add large particles
        for(auto iterLarger = _hierarchyLevels[largerLevel].begin(); iterLarger.isValid(); ++iterLarger) {
          const double particleRad = - 1.0 * iterLarger->getRad();
          iterLarger->setRad(particleRad);
          crossLevel.addParticleImpl(*iterLarger);
          }
        _hierarchyLevels[largerLevel].deleteAllParticles();
        

        //Add small particles
        for(auto iterSmaller = _hierarchyLevels[smallerLevel].begin(); iterSmaller.isValid(); ++iterSmaller) {
          crossLevel.addParticleImpl(*iterSmaller);
        }
        _hierarchyLevels[smallerLevel].deleteAllParticles();

        // Iterate pairwise over the cross-level
        crossLevel.iteratePairwise(traversal);

        // Move particles back to their level
        AUTOPAS_OPENMP(parallel)
        for (auto iterCross = crossLevel.begin(); iterCross.isValid(); ++iterCross) {

          const double particleRad = iterCross->getRad();

          //if targetLevel = largerLevel, particle has negative radius. Hence, reverse the ugly trick
          if (particleRad < 0.) {
            
            iterCross->setRad(std::abs(particleRad));
          }

          unsigned int targetLevel = getHierarchyLevelOfParticle(*iterCross);

          addParticleToGridLevel(*iterCross, targetLevel);
        }
        
        crossLevel.deleteAllParticles();

      }
    }

    demLib::DEMFunctor<Particle>::endCrossLevelTraversal();
  }

  /**
   * @todo
   * 
   * Iterates over all particle pairs in the container.
   * @param traversal The traversal to use for the iteration.
   */
  void iteratePairwise(TraversalInterface *traversal) override {

    // Distribute particles over the hierarchical grid levels, if more than one HG level is present
    if (_hierarchyLevels.size() > 1) {
      
      distributeParticlesToHGLevels();

    }
    // Else skip the distribution and just iteratePairwise over the carrier level and return
    else {
      
      _hierarchyLevels[0].iteratePairwise(traversal);

      return;

    }

    // Iterate over hierarchy levels 
    iterateAllHGLevels(traversal);

    // Iterate over cross levels 
    iterateAllCrossLevels(traversal);

    
    // Distribute particles back to the top level hierarchical grid
    distributeParticlesToCarrierLevel();

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
   * @copydoc autopas::ParticleContainerInterface::updateContainer()
   */
  [[nodiscard]] std::vector<ParticleType> updateContainer(bool keepNeighborListsValid) override {

    return _hierarchyLevels[0].updateContainer(keepNeighborListsValid);

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {

    return _hierarchyLevels[0].getTraversalSelectorInfo();

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getParticle()
   */
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                            IteratorBehavior iteratorBehavior) const override {

    return _hierarchyLevels[0].getParticle(cellIndex, particleIndex, iteratorBehavior); 

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getParticle()
   */
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                                   IteratorBehavior iteratorBehavior,
                                                                   const std::array<double, 3> &boxMin,
                                                                   const std::array<double, 3> &boxMax) const override {

    return _hierarchyLevels[0].getParticle(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::deleteParticle()
   */
  bool deleteParticle(Particle &particle) override {

    return _hierarchyLevels[0].deleteParticle(particle);

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::deleteParticle()
   */
  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {

    return _hierarchyLevels[0].deleteParticle(cellIndex, particleIndex);

  }

  /**
   * @brief Get the Hierarchy Level Of Particle object
   * 
   * @param p 
   * @return Level the particle should be sorted into. 
   */
  unsigned int getHierarchyLevelOfParticle(const ParticleType &p) {

    // If only one HG level is available all particles are sorted into this level.
    if (_hierarchyLevels.size() == 0) {
      return 0;
    }
        
    // Get the radius of the particle in question
    auto particleRad = p.getRad();

    unsigned int level;

    // Loop over all hierarchy levels
    for (level = 0; level < _hierarchyLevels.size(); level++) {
      // As long as the '_hierarchicalGridBorders' element is smaller than the radius continue to cycle
      if (particleRad > _hierarchicalGridBorders[level]) break;     
    }

    if(level == 0) {
    utils::ExceptionHandler::exception(
      "HierarchicalGrids: Trying to get the hierarchy level of a particle which is too large to fit in the grid.\n"
      "Particle radius {}\n"
      "Upper limit of the largest level {}",
      particleRad, _hierarchicalGridBorders[level]);
    }
    // Break and return the level
    return level - 1;
  }

  void setHierarchicalGridBorders(std::vector<double> boundaries) {

    _hierarchicalGridBorders = boundaries;

  }

  std::vector<double> getHierarchicalGridBorders() {

    return _hierarchicalGridBorders;

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
  * @note This can not be a std::vector<LinkedCells<Particle>> because there is no copy constructor for LinkedCells. 
  */
  std::deque<LinkedCells<Particle>> _hierarchyLevels;

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
   * By design '_hierarchicalGridBorders' should have the same number of levels as the H-Grid.
   * 
   */
  std::vector<double> _hierarchicalGridBorders;

};
  
} // namespace autopas 
