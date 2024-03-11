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
  using ParticleCell = FullParticleCell<Particle>;

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

    _crossLevelInteractions.reserve(getNumberOfCrossLevelInteractions());

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

    // Get level of the particle and add it to the corresponding level
    unsigned int level = getHierarchyLevelOfParticle(p);
    addParticleToGridLevel(p, level);

  }

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

  /**
   * @brief Adds a particle to the container that lies in the halo region of the container.
   * 
   * This is an unsafe version of addHaloParticle() and does not perform a boundary check.
   * @param haloParticle Particle to be added. This particle is already checked to be outside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be outside of the bounding box!
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {

    // Get level of halo particle and add it to the corresponding level
    unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    _hierarchyLevels[level].addHaloParticle(haloParticle);

  }

  public:
  /** 
   * @copydoc autopas::ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {

    // Get the level of the halo particle and update
    unsigned int level = getHierarchyLevelOfParticle(haloParticle);
    return _hierarchyLevels[level].updateHaloParticle(haloParticle);

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
    for (auto & hgLevel : _hierarchyLevels) {
      hgLevel.deleteHaloParticles();
    }
  }

  /**
   * Deletes all particles.
   */
  void deleteAllParticles() override {
    for (auto & hgLevel : _hierarchyLevels) {
      hgLevel.deleteAllParticles();
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getNumberOfParticles()
   */
  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const override {
    size_t numParticles = 0ul;
    for (auto & hgLevel : _hierarchyLevels) {
      numParticles += hgLevel.getNumberOfParticles(behavior);
    }
    return numParticles;
    
  }

  /**
   *@copydoc autopas::ParticleContainerInterface::size()
   */
  [[nodiscard]] size_t size() const override {
    size_t numParticles = 0ul;
    for (auto & hgLevel : _hierarchyLevels) {
      numParticles += hgLevel.size();
    }
    return numParticles;
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<ParticleType, true, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, true, false>::ParticleVecType *additionalVectors = nullptr) override {
        
    return ContainerIterator<ParticleType, true, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note const version
   */
  [[nodiscard]] ContainerIterator<ParticleType, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, false, false>::ParticleVecType *additionalVectors = nullptr) const override {
        
    return ContainerIterator<ParticleType, false, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<ParticleType, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, true, true>::ParticleVecType *additionalVectors = nullptr) override {

    return ContainerIterator<ParticleType, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   * @note const version
   */
  [[nodiscard]] ContainerIterator<ParticleType, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, false, true>::ParticleVecType *additionalVectors = nullptr) const override {

    return ContainerIterator<ParticleType, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc autopas::LinkedCells::forEachInRegion()
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    
    // Apply Lambda to each level 
    for (auto & hgLevel : _hierarchyLevels) {
      hgLevel.forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    }

  }

  /**
   * @copydoc autopas::LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    
    // Apply reduction to each level
    for (auto & hgLevel : _hierarchyLevels) {
      hgLevel.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    }

  } 

  protected:

  void iterateAllHGLevels(TraversalInterface *traversal) {

    // @idea: Reduce the number of threads per level and parallelize as much as 
    // possible over the levels. Might mean some major redesign all over autopas? 
    //AUTOPAS_OPENMP(parallel firstprivate (traversal))
    for (auto & hgLevel : _hierarchyLevels) {

      hgLevel.iteratePairwise(traversal);
      
    }
  }

  void iterateAllCrossLevels() {

    AUTOPAS_OPENMP(parallel)
    for (auto & levelPairs : _crossLevelInteractions) {

      size_t largerLevel, smallerLevel;
      std::tie(largerLevel, smallerLevel) = levelPairs;

      iterateOverCrossLevels(largerLevel, smallerLevel);

    }
    
  }

  void iterateOverCrossLevels(unsigned int firstLevel, unsigned int secondLevel) {


    demLib::DEMFunctor<Particle> demFunctor(1.);

    //iterate over all particles in a given larger level
    // @note Maybe parallelize iteration over second level?
    //AUTOPAS_OPENMP(parallel firstprivate(demFunctor)) // <- cannot parallelize here, since this leads to erroneous forces
    for(auto iterFirstLVLParticles = _hierarchyLevels[firstLevel].begin(); iterFirstLVLParticles.isValid(); ++iterFirstLVLParticles) {

      // Get position and radius of larger particle
      const std::array<double, 3> posLarger = iterFirstLVLParticles->getR();
      const double radLarger = iterFirstLVLParticles->getRad();

      // Box size of particle interations
      const double crossLevelCellSize = ( radLarger + _hierarchicalGridBorders[secondLevel] );

      // Get lower and higher corner of box around the larger particle containing possible interaction partners
      // in the smaller level
      const std::array<double, 3> lowerBoxCorner = autopas::utils::ArrayMath::subScalar(posLarger, crossLevelCellSize);
      const std::array<double, 3> higherBoxCorner = autopas::utils::ArrayMath::addScalar(posLarger, crossLevelCellSize);
      
      AUTOPAS_OPENMP(parallel firstprivate(demFunctor))
      for (auto iterSecondLVLParticles = 
          _hierarchyLevels[secondLevel].getRegionIterator(lowerBoxCorner, higherBoxCorner, autopas::IteratorBehavior::ownedOrHalo); 
                iterSecondLVLParticles.isValid(); ++iterSecondLVLParticles) {

        // Calculate forces
        // @note Maybe use different functor?
        demFunctor.AoSFunctor(*iterFirstLVLParticles, *iterSecondLVLParticles, true);
      }

    }

  }

  public:

  /**
   * @todo
   * 
   * Iterates over all particle pairs in the container.
   * @param traversal The traversal to use for the iteration.
   */
  void iteratePairwise(TraversalInterface *traversal) override {

    // Iterate over hierarchy levels 
    iterateAllHGLevels(traversal);

    // Iterate over cross levels 
    // @note No traversal needed, since forces between HG-levels are calculated directly
    iterateAllCrossLevels();

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

    // Set up vector collecting all particles returned by applying updateContainer()
    // to the HG-levels
    std::vector<ParticleType> collectParticles;
    // TODO: needs smarter heuristic than this. Taken from autopas::LinkedCells::updateContainer()
    collectParticles.reserve(128);

    for (auto & hgLevel : _hierarchyLevels) {
      // Iterate over the HG-levels, update each level and collect the returned particles
      std::vector<ParticleType> collectParticlesOfLevel = hgLevel.updateContainer(keepNeighborListsValid);

      // Merge with global vector
      collectParticles.insert(collectParticles.end(), collectParticlesOfLevel.begin(),
                                                      collectParticlesOfLevel.end());
    }
    

    return collectParticles;

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
                                                                   IteratorBehavior iteratorBehavior,
                                                                   const std::array<double, 3> &boxMin,
                                                                   const std::array<double, 3> &boxMax) const override {

    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getParticle()
   */
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                            IteratorBehavior iteratorBehavior) const override {

    // this is not a region iter hence we stretch the bounding box to the numeric max
    constexpr std::array<double, 3> boxMin{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                                           std::numeric_limits<double>::lowest()};

    constexpr std::array<double, 3> boxMax{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                           std::numeric_limits<double>::max()};
    return getParticleImpl<false>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
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
  std::tuple<const Particle *, size_t, size_t> getParticleImpl(size_t encodedLevelAndCellIndex, size_t particleIndex,
                                                               IteratorBehavior iteratorBehavior,
                                                               const std::array<double, 3> &boxMin,
                                                               const std::array<double, 3> &boxMax) const {
    
    const Particle *retPtr = nullptr;
    size_t currentCellIndex;
    size_t currentParticleIndex;
    size_t currentLevelAndCellIndex;

    //decode HG-level and cellIndex
    auto [currentLevel, cellIndex] = decodeCellAndLevel(encodedLevelAndCellIndex);

    //and get the particle in the corresponding HG-level 
    if constexpr (regionIter) {
      std::tie(retPtr, currentCellIndex, currentParticleIndex) = 
          _hierarchyLevels[currentLevel].getParticle(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
    } else {
      std::tie(retPtr, currentCellIndex, currentParticleIndex) = 
          _hierarchyLevels[currentLevel].getParticle(cellIndex, particleIndex, iteratorBehavior);
    }

    // Check where we are
    // If _hierarchyLevels[currentLevel].getParticle(...) returned a nullptr, we are out of
    // bounds on the current level of the hierarchical grid.
    if (retPtr == nullptr) { 
      // If the current level was the smallest grid level, we are done iterating
      if (currentLevel == _hierarchyLevels.size() - 1) {

        return {nullptr, 0, 0};

      } else {
      // If currentLevel was not the smallest grid level, we finished iterating over this
      // specific level and we have to advance to following one.
      
      // Increase current level by one
      currentLevel += 1;

      // Obtain the first particle of the following level
      if constexpr (regionIter) {
        std::tie(retPtr, currentCellIndex, currentParticleIndex) = 
            _hierarchyLevels[currentLevel].getParticle(0, 0, iteratorBehavior, boxMin, boxMax);
      } else {
        std::tie(retPtr, currentCellIndex, currentParticleIndex) = 
            _hierarchyLevels[currentLevel].getParticle(0, 0, iteratorBehavior);
      }

      // Encode the level and cell index
      currentLevelAndCellIndex = encodeCellAndLevel(currentLevel, currentCellIndex);

      // And return for the following iteration
      return {retPtr, currentLevelAndCellIndex, currentParticleIndex};

      }
    }

    // If retPtr is NOT a nullptr, we are still iterating over the current level 
    //
    // Encode the level and cell index
    currentLevelAndCellIndex = encodeCellAndLevel(currentLevel, currentCellIndex);

    // And return for the following iteration
    return {retPtr, currentLevelAndCellIndex, currentParticleIndex};
  }
  
  /**
   * @copydoc autopas::ParticleContainerInterface::deleteParticle()
   */
  bool deleteParticle(Particle &particle) override {

    // get level of particle and delete it from there
    size_t level = getHierarchyLevelOfParticle(particle);
    return _hierarchyLevels[level].deleteParticle(particle);

  }

  /**
   * @copydoc autopas::ParticleContainerInterface::deleteParticle()
   * 
   * @note This version only works correctly if the cellIndex is actually an 
   * encoded cell and level index
   */
  bool deleteParticle(size_t encodedCellAndLevelIndex, size_t particleIndex) override {

    auto [level, cellIndex] = decodeCellAndLevel(encodedCellAndLevelIndex);

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
   * Encodes a level number and a cell index as a combined value. The logic behind this is
   * encodedCellAndLevel = cellIndex + _levelOffset * level
   * with _levelOffset = 1e15 (see below).
   * 
   * @note This fails if the number of cells in a level is larger than _levelOffset
   * 
   * @param level 
   * @param cellIndex 
   * @return size_t Encoded cell and level
   */
  size_t encodeCellAndLevel(size_t level, size_t cellIndex) const {

    // Exception if current level contains too many cells
    if (cellIndex >= _levelOffset) {
      utils::ExceptionHandler::exception(
        "HierarchicalGrids: Trying to encode a cellIndex which is larger than _levelOffset = 1e15,\n"
        "i.e. the current HG-level contains too many cells.\n"
        "Current level: {}", level);
    }
    
    
    return cellIndex + _levelOffset * level;
  }

  /**
   * Decode a encoded cellID and HG level.
   * 
   * @param encodedCellPlusLevel 
   * @return std::tuple<size_t, size_t> level, cellID
   */
  std::tuple<size_t, size_t> decodeCellAndLevel(size_t encodedCellPlusLevel) const {

    size_t level = encodedCellPlusLevel / _levelOffset;
    size_t cellID = encodedCellPlusLevel - level * _levelOffset;

    return {level, cellID};
  }

  /**
   * Calculate the number of cross-level interactions for the present number of hierarchy levels
   * 
   * @return size_t number of cross-level interactions
   */
  size_t getNumberOfCrossLevelInteractions() const {

    size_t numberOfCrossLevels = 0;

    for (unsigned int largerLevel = 0; largerLevel < _hierarchyLevels.size(); largerLevel++) {
      for (unsigned int smallerLevel = largerLevel + 1; smallerLevel < _hierarchyLevels.size(); smallerLevel++) {

        numberOfCrossLevels++;

      }
    }

    return numberOfCrossLevels;

  }

  /**
   * Set the cross-level interaction pairs, i.e. the which larger level is checked with which smaller level
   * 
   */
  void setCrossLevelInteractionPairs() {

    for (unsigned int largerLevel = 0; largerLevel < _hierarchyLevels.size(); largerLevel++) {
      for (unsigned int smallerLevel = largerLevel + 1; smallerLevel < _hierarchyLevels.size(); smallerLevel++) {

        _crossLevelInteractions.emplace_back(largerLevel, smallerLevel);

      }
    }

  }



  /**
  * @brief Vector containing the HGrid's different hierarchy levels of Linked Cells.
  * 
  * @note This can not be a std::vector<LinkedCells<Particle>> because there is no copy constructor for LinkedCells. 
  */
  std::deque<LinkedCells<Particle>> _hierarchyLevels;

  std::vector<std::tuple<size_t, size_t>> _crossLevelInteractions;

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

  /**
   * Level offset required for encoding and decoding the current HG-level and the cellIndex 
   * @see HierarchicalGrids::getParticleImpl.
   * 
   */
  const size_t _levelOffset = 1e15;

};
  
} // namespace autopas 
