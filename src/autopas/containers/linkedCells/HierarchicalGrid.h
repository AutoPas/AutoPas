/**
 * @file HierarchicalGrid.h
 *
 * @date 29.11.2024
 * @author atacann
 */
#pragma once

#include <algorithm>
#include <array>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/cellTraversals/BalancedTraversal.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * HierarchicalGrid class.
 * This class stores multiple LinkedCells containers, each for different sizes of particles
 * Traverses them all one by one, and cross hierarchy afterward
 * @tparam Particle type of the Particle
 */
template <class Particle>
class HierarchicalGrid : public ParticleContainerInterface<Particle> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle;

    /**
   * Constructor of the HierarchicalGrid class.
   * @param boxMin
   * @param boxMax
   * @param baseCutoff base cutoff for each particle, it will be scaled by the size of a Particle
   * @param cutoffs cutoffs for each level of the hierarchy
   * @param skinPerTimestep
   * @param rebuildFrequency
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param loadEstimator the load estimation algorithm for balanced traversals.
   */
  HierarchicalGrid(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double baseCutoff, 
              const std::vector<double>& cutoffs, const double skinPerTimestep, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0)
      : ParticleContainerInterface<Particle>(skinPerTimestep),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _baseCutoff(baseCutoff),
        _skin(skinPerTimestep * rebuildFrequency),
        _numHierarchyLevels(cutoffs.size()),
        _cellSizeFactor(cellSizeFactor),
        _cutoffs(cutoffs) {
          if(_cutoffs.empty()) {
            utils::ExceptionHandler::exception("Error: Hierarchical Grid cutoffs vector is empty.");
          }
          // make sure cutoffs are sorted
          std::sort(_cutoffs.begin(), _cutoffs.end());
          // biggest interaction length
          double interactionLengthBiggest = _cutoffs.back() + _skin;
          // generate LinkedCells for each hierarchy, with different cutoffs
          _hierarchies.reserve(_numHierarchyLevels);
          for(size_t i = 0; i < _numHierarchyLevels; ++i) {
            // here, LinkedCells with smaller cutoffs need to have bigger halo regions, because particles on level above
            // can potentially interact with those, even though these halo particles cannot interact with other particles inside its own level
            // a hacky way: make all LinkedCells interactionLength equal to the biggest one, adjust cellSizeFactor
            // for each LinkedCells so that the cellLength is equal to actual interaction length of that level
            const double interactionLengthLevel = _cutoffs[i] + _skin;
            const double ratio = interactionLengthLevel / interactionLengthBiggest;
            _hierarchies[i] = LinkedCells<Particle>(_boxMin, _boxMax, _cutoffs.back(), skinPerTimestep, rebuildFrequency, cellSizeFactor * ratio, loadEstimator);
          }
        }

  /**
 * Destructor of HierarchicalGrid.
 */
  ~HierarchicalGrid() override = default;

  /**
   * Delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  HierarchicalGrid(const HierarchicalGrid &obj) = delete;

  /**
   * Delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  HierarchicalGrid &operator=(const HierarchicalGrid &other) = delete;

  /**
   * @copydoc autopas::ParticleContainerInterface::getCutoff()
   */
  [[nodiscard]] double getCutoff() const final { return _baseCutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) final { _baseCutoff = cutoff; }

  [[nodiscard]] double getVerletSkin() const final { return _skin; }

  [[nodiscard]] virtual size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const override {
    size_t numParticles = 0;
    for (const auto &linkedCells: _hierarchies) {
      numParticles += linkedCells.getNumberOfParticles();
    }
    return numParticles;
  }

  [[nodiscard]] virtual size_t size() const override {
    size_t numParticles = 0;
    for (const auto &linkedCells: _hierarchies) {
      numParticles += linkedCells.size();
    }
    return numParticles;
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
 * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
 * In HierarchicalGrids, interaction length returns maximum interaction length which is the biggest cutoff + skin
 */
  [[nodiscard]] double getInteractionLength() const final {
    return _cutoffs.back() + _skin;
  }
  
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::hierarchicalGrid; }

  [[nodiscard]] CellType getParticleCellTypeEnum() const override { return CellType::FullParticleCell; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    // TODO: add reserve somehow
    // do not do anything, because numParticles includes all particles over all hierarchies
    // can do numParticles / _numHierarchyLevels, will not be accurate
  }

  void addParticleImpl(const ParticleType &p) override {
    _hierarchies[getHierarchyLevel(p)].addParticle(p);
  }

  bool deleteParticle(Particle &particle) override {
    return _hierarchies[getHierarchyLevel(particle)].deleteParticle(particle);
  }

  // TODO: Implement
  virtual bool deleteParticle(size_t cellIndex, size_t particleIndex) override = 0;

  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    _hierarchies[getHierarchyLevel(haloParticle)].addParticle(haloParticle);
  }

  bool updateHaloParticle(const ParticleType &haloParticle) override {
    return _hierarchies[getHierarchyLevel(haloParticle)].updateHaloParticle(haloParticle);
  }

  void deleteHaloParticles() override {
    AUTOPAS_OPENMP(parallel for)
    for (size_t idx = 0; idx < _hierarchies.size(); idx++) {
      _hierarchies[idx].deleteHaloParticles();
    }
  }

  void deleteAllParticles() override {
    for (auto &linkedCells: _hierarchies) {
      linkedCells.deleteAllParticles();
    }
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
  }

  void computeInteractions(TraversalInterface *traversal) override {
    // calculate only within individual LinkedCells
    // I think not parallelizing here is better, it is already parallelized in computeInteractions
    for (size_t i = 0; i < _hierarchies.size(); i++) {
      // here traversal needs to be built again with the correct interactionLength, as it will be wrong (except for the highest level)
      const auto traversalOption = traversal->getTraversalType();
      const auto newton3 = traversal->getUseNewton3();
      const auto dataLayout = traversal->getDataLayout();
      const auto traversalSelectorInfo = _hierarchies[i].getTraversalSelectorInfo();
      // fix interactionLength
      traversalSelectorInfo.interactionLength = _cutoffs[i] + _skin; // cellLength
      // TODO: how to get functor here?
      _hierarchies[i].computeInteractions(traversal);
    }
    // calculate cross hierarchy
    computeCrossInteractions(traversal);
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return autopas::LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }
    std::vector<ParticleType> invalidParticles;
    // parallelized inside each updateContainer
    for (auto &linkedCells: _hierarchies) {
      auto invalidParticlesSingle = linkedCells.updateContainer(keepNeighborListsValid);
      invalidParticles.insert(invalidParticles.end(), std::make_move_iterator(invalidParticlesSingle.begin()),
        std::make_move_iterator(invalidParticlesSingle.end()));
    }
    return invalidParticles;
  }

  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // return traversal info of hierarchy with biggest interactionLength
    return _hierarchies.back().getTraversalSelectorInfo();
  }

  // TODO: Implement all below
 //  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
 //                                                         IteratorBehavior iteratorBehavior,
 //                                                         const std::array<double, 3> &boxMin,
 //                                                         const std::array<double, 3> &boxMax) const override {
 //    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
 //  }
 //  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
 //                                                           IteratorBehavior iteratorBehavior) const override {
 //    // this is not a region iter hence we stretch the bounding box to the numeric max
 //    constexpr std::array<double, 3> boxMin{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
 //                                           std::numeric_limits<double>::lowest()};
 //
 //    constexpr std::array<double, 3> boxMax{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
 //                                           std::numeric_limits<double>::max()};
 //    return getParticleImpl<false>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
 //  }
 //
 //    /**
 //   * Container specific implementation for getParticle. See ParticleContainerInterface::getParticle().
 //   *
 //   * @tparam regionIter
 //   * @param cellIndex
 //   * @param particleIndex
 //   * @param iteratorBehavior
 //   * @param boxMin
 //   * @param boxMax
 //   * @return tuple<ParticlePointer, CellIndex, ParticleIndex>
 //   */
 //  template <bool regionIter>
 //  std::tuple<const Particle *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
 //                                                               IteratorBehavior iteratorBehavior,
 //                                                               const std::array<double, 3> &boxMin,
 //                                                               const std::array<double, 3> &boxMax) const {
 //    using namespace autopas::utils::ArrayMath::literals;
 //
 //    std::array<double, 3> boxMinWithSafetyMargin = boxMin;
 //    std::array<double, 3> boxMaxWithSafetyMargin = boxMax;
 //    if constexpr (regionIter) {
 //      // We extend the search box for cells here since particles might have moved
 //      boxMinWithSafetyMargin -= (this->_skinPerTimestep * static_cast<double>(this->getStepsSinceLastRebuild()));
 //      boxMaxWithSafetyMargin += (this->_skinPerTimestep * static_cast<double>(this->getStepsSinceLastRebuild()));
 //    }
 //
 //    // first and last relevant cell index
 //    const auto [startCellIndex, endCellIndex] = [&]() -> std::tuple<size_t, size_t> {
 //      if constexpr (regionIter) {
 //        // We extend the search box for cells here since particles might have moved
 //        return {_cellBlock.get1DIndexOfPosition(boxMinWithSafetyMargin),
 //                _cellBlock.get1DIndexOfPosition(boxMaxWithSafetyMargin)};
 //      } else {
 //        if (not(iteratorBehavior & IteratorBehavior::halo)) {
 //          // only potentially owned region
 //          return {_cellBlock.getFirstOwnedCellIndex(), _cellBlock.getLastOwnedCellIndex()};
 //        } else {
 //          // whole range of cells
 //          return {0, this->_cells.size() - 1};
 //        }
 //      }
 //    }();
 //
 //    // if we are at the start of an iteration ...
 //    if (cellIndex == 0 and particleIndex == 0) {
 //      cellIndex =
 //          startCellIndex + ((iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num());
 //    }
 //    // abort if the start index is already out of bounds
 //    if (cellIndex >= this->_cells.size()) {
 //      return {nullptr, 0, 0};
 //    }
 //    // check the data behind the indices
 //    if (particleIndex >= this->_cells[cellIndex].size() or
 //        not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
 //            this->_cells[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax)) {
 //      // either advance them to something interesting or invalidate them.
 //      std::tie(cellIndex, particleIndex) =
 //          advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax,
 //                                             boxMinWithSafetyMargin, boxMaxWithSafetyMargin, endCellIndex);
 //    }
 //
 //    // shortcut if the given index doesn't exist
 //    if (cellIndex > endCellIndex) {
 //      return {nullptr, 0, 0};
 //    }
 //    const Particle *retPtr = &this->_cells[cellIndex][particleIndex];
 //
 //    return {retPtr, cellIndex, particleIndex};
 //  }
 //
 //  [[nodiscard]] ContainerIterator<ParticleType, true, false> begin(
 //    IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
 //    typename ContainerIterator<ParticleType, true, false>::ParticleVecType *additionalVectors = nullptr) override {
 //    return ContainerIterator<ParticleType, true, false>(*this, behavior, additionalVectors);
 //  }
 //
 //  [[nodiscard]] ContainerIterator<ParticleType, false, false> begin(
 //      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
 //      typename ContainerIterator<ParticleType, false, false>::ParticleVecType *additionalVectors =
 //          nullptr) const override {
 //    return ContainerIterator<ParticleType, false, false>(*this, behavior, additionalVectors);
 //  }
 //
 //    /**
 //   * Execute code on all particles in this container in a certain region as defined by a lambda function.
 //   * @tparam Lambda (Particle &p) -> void
 //   * @param forEachLambda code to be executed on all particles
 //   * @param lowerCorner lower corner of bounding box
 //   * @param higherCorner higher corner of bounding box
 //   * @param behavior @see IteratorBehavior
 //   */
 //  template <typename Lambda>
 //  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
 //                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
 //    using namespace autopas::utils::ArrayMath::literals;
 //
 //    const auto startIndex3D = this->_cellBlock.get3DIndexOfPosition(lowerCorner - this->getVerletSkin());
 //    const auto stopIndex3D = this->_cellBlock.get3DIndexOfPosition(higherCorner + this->getVerletSkin());
 //
 //    const size_t numCellsOfInterest = (stopIndex3D[0] - startIndex3D[0] + 1) * (stopIndex3D[1] - startIndex3D[1] + 1) *
 //                                      (stopIndex3D[2] - startIndex3D[2] + 1);
 //    std::vector<size_t> cellsOfInterest;
 //    cellsOfInterest.reserve(numCellsOfInterest);
 //
 //    const auto &cellsPerDimensionWithHalo = this->_cellBlock.getCellsPerDimensionWithHalo();
 //
 //    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
 //      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
 //        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
 //          cellsOfInterest.push_back(utils::ThreeDimensionalMapping::threeToOneD({x, y, z}, cellsPerDimensionWithHalo));
 //        }
 //      }
 //    }
 //
 //    for (auto cellIndex : cellsOfInterest) {
 //      if (not _cellBlock.ignoreCellForIteration(cellIndex, behavior)) {
 //        getCells()[cellIndex].forEach(forEachLambda, lowerCorner, higherCorner, behavior);
 //      }
 //    }
 //  }
 //
 //  /**
 //   * Execute code on all particles in this container in a certain region as defined by a lambda function.
 //   * @tparam Lambda (Particle &p, A &result) -> void
 //   * @tparam A type of reduction Value
 //   * @param reduceLambda code to be executed on all particles
 //   * @param result reference to starting and final value for reduction
 //   * @param lowerCorner lower corner of bounding box
 //   * @param higherCorner higher corner of bounding box
 //   * @param behavior @see IteratorBehavior
 //   */
 //  template <typename Lambda, typename A>
 //  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
 //                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
 //    using namespace autopas::utils::ArrayMath::literals;
 //    const auto startIndex3D = this->_cellBlock.get3DIndexOfPosition(lowerCorner - this->getVerletSkin());
 //    const auto stopIndex3D = this->_cellBlock.get3DIndexOfPosition(higherCorner + this->getVerletSkin());
 //
 //    const size_t numCellsOfInterest = (stopIndex3D[0] - startIndex3D[0] + 1) * (stopIndex3D[1] - startIndex3D[1] + 1) *
 //                                      (stopIndex3D[2] - startIndex3D[2] + 1);
 //    std::vector<size_t> cellsOfInterest;
 //    cellsOfInterest.reserve(numCellsOfInterest);
 //
 //    const auto &cellsPerDimensionWithHalo = this->_cellBlock.getCellsPerDimensionWithHalo();
 //
 //    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
 //      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
 //        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
 //          cellsOfInterest.push_back(utils::ThreeDimensionalMapping::threeToOneD({x, y, z}, cellsPerDimensionWithHalo));
 //        }
 //      }
 //    }
 //
 //    for (auto cellIndex : cellsOfInterest) {
 //      if (not _cellBlock.ignoreCellForIteration(cellIndex, behavior)) {
 //        getCells()[cellIndex].reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
 //      }
 //    }
 //  }
 //
 //  /**
 // * Get the cell block, not supposed to be used except by verlet lists
 // * @return the cell block
 // */
 //  internal::CellBlock3D<ParticleCell> &getCellBlock() { return nullptr; }
 //
 //  /**
 // * @copydoc autopas::LinkedCells::getCellBlock()
 // * @note const version
 // */
 //  const internal::CellBlock3D<ParticleCell> &getCellBlock() const { return nullptr; }
 //
 //  /**
 //   * Returns a non-const reference to the cell data structure.
 //   * @return Non-const reference.
 //   */
 //  std::vector<ParticleCell> &getCells() { return nullptr; }
 //
 //  /**
 //  * Get immutable vector of cells.
 //  * @return immutable reference to _cells
 //  */
 //  [[nodiscard]] const std::vector<ParticleCell> &getCells() const {
 //    return std::vector<ParticleCell>();
 //  }


 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _baseCutoff;
  double _skin;
  size_t _numHierarchyLevels;
  const double _cellSizeFactor;

  /**
   *
   * @param p Particle to add into HierarchicalGrid
   * @return which Hierarchy the particle belongs to
   */
  size_t getHierarchyLevel(const ParticleType &p) const {
    // binary search not worth if there are small amount of levels
    // scale size by baseCutoff
    const double cutoff = p.getSize() * _baseCutoff;
    for (size_t i = 0; i < _numHierarchyLevels; ++i) {
      if (_cutoffs[i] >= cutoff) {
        return i;
      }
    }
    AutoPasLog(ERROR, "Size of Particle times baseCutoff is bigger than biggest cutoff of HierarchicalGrid, "
                      "will result in wrong interaction calculation");
    return _numHierarchyLevels - 1;
  }

  /**
   * Computer interactions across hierachies.
   * @param traversal Not used
   */
  void computeCrossInteractions(TraversalInterface *traversal) {
    // TODO: implement
  }

 protected:
  //   /**
  //  * Given a pair of cell-/particleIndex and iterator restrictions either returns the next indices that match these
  //  * restrictions or indices that are out of bounds (e.g. cellIndex >= cells.size())
  //  * @tparam regionIter
  //  * @param cellIndex
  //  * @param particleIndex
  //  * @param iteratorBehavior
  //  * @param boxMin The actual search box min
  //  * @param boxMax The actual search box max
  //  * @param boxMinWithSafetyMargin Search box min that includes a surrounding of skinPerTimestep * stepsSinceLastRebuild
  //  * @param boxMaxWithSafetyMargin Search box max that includes a surrounding of skinPerTimestep * stepsSinceLastRebuild
  //  * @param endCellIndex Last relevant cell index
  //  * @return tuple<cellIndex, particleIndex>
  //  */
  // template <bool regionIter>
  // std::tuple<size_t, size_t> advanceIteratorIndices(
  //     size_t cellIndex, size_t particleIndex, IteratorBehavior iteratorBehavior, const std::array<double, 3> &boxMin,
  //     const std::array<double, 3> &boxMax, std::array<double, 3> boxMinWithSafetyMargin,
  //     std::array<double, 3> boxMaxWithSafetyMargin, size_t endCellIndex) const {
  //   // Finding the indices for the next particle
  //   const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();
  //
  //   // helper function to determine if the cell can even contain particles of interest to the iterator
  //   auto cellIsRelevant = [&]() -> bool {
  //     bool isRelevant =
  //         // behavior matches possible particle ownership
  //         (iteratorBehavior & IteratorBehavior::owned and _cellBlock.cellCanContainOwnedParticles(cellIndex)) or
  //         (iteratorBehavior & IteratorBehavior::halo and _cellBlock.cellCanContainHaloParticles(cellIndex));
  //     if constexpr (regionIter) {
  //       // short circuit if already false
  //       if (isRelevant) {
  //         // is the cell in the region?
  //         const auto [cellLowCorner, cellHighCorner] = _cellBlock.getCellBoundingBox(cellIndex);
  //         isRelevant =
  //             utils::boxesOverlap(cellLowCorner, cellHighCorner, boxMinWithSafetyMargin, boxMaxWithSafetyMargin);
  //       }
  //     }
  //     return isRelevant;
  //   };
  //
  //   do {
  //     // advance to the next particle
  //     ++particleIndex;
  //     // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.
  //
  //     // If cell has wrong type, or there are no more particles in this cell jump to the next
  //     while (not cellIsRelevant() or particleIndex >= this->_cells[cellIndex].size()) {
  //       // TODO: can this jump be done more efficient if behavior is only halo or owned?
  //       // TODO: can this jump be done more efficient for region iters if the cell is outside the region?
  //       cellIndex += stride;
  //       particleIndex = 0;
  //
  //       // If we notice that there is nothing else to look at set invalid values, so we get a nullptr next time and
  //       // break.
  //       if (cellIndex > endCellIndex) {
  //         return {std::numeric_limits<decltype(cellIndex)>::max(), std::numeric_limits<decltype(particleIndex)>::max()};
  //       }
  //     }
  //   } while (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
  //       this->_cells[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax));
  //
  //   // the indices returned at this point should always be valid
  //   return {cellIndex, particleIndex};
  // }
  //
  // /**
  //  * Checks if a given traversal is allowed for LinkedCells and sets it up for the force interactions.
  //  * @tparam Traversal Traversal type. E.g. pairwise, triwise
  //  * @param traversal
  //  */
  // template <typename Traversal>
  // void prepareTraversal(Traversal &traversal) {
  //   auto *traversalInterface = dynamic_cast<LCTraversalInterface *>(traversal);
  //   auto *cellTraversal = dynamic_cast<CellTraversal<ParticleCell> *>(traversal);
  //   // if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
  //   //   balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
  //   // }
  //   if (traversalInterface && cellTraversal) {
  //     cellTraversal->setCellsToTraverse(this->_cells);
  //   } else {
  //     autopas::utils::ExceptionHandler::exception(
  //         "The selected traversal is not compatible with the LinkedCells container. TraversalID: {}",
  //         traversal->getTraversalType());
  //   }
  // }

  std::vector<autopas::LinkedCells<Particle>> _hierarchies;
  std::vector<double> _cutoffs;

};

}