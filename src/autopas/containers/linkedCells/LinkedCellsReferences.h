/**
 * @file LinkedCells.h
 * @date 05.04.2020
 * @author lunaticcoding
 */

#pragma once

#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/containers/cellTraversals/BalancedTraversal.h"
#include "autopas/containers/linkedCells/ParticleVector.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ParticleCellHelpers.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "autopas/utils/markParticleAsDeleted.h"

namespace autopas {

/**
 * LinkedCells class.
 * This class uses a list of neighboring cells to store the particles.
 * These cells dimensions are at least as large as the given cutoff radius,
 * therefore short-range interactions only need to be calculated between
 * particles in neighboring cells.
 * @tparam Particle_T type of the Particle
 */
template <class Particle_T>
class LinkedCellsReferences : public CellBasedParticleContainer<ReferenceParticleCell<Particle_T>> {
 public:
  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle_T;
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCellType = ReferenceParticleCell<Particle_T>;

  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param sortingThreshold number of particles in two cells from which sorting should be performed
   * @param loadEstimator the load estimation algorithm for balanced traversals.
   * By default all applicable traversals are allowed.
   */
  LinkedCellsReferences(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
                        const double skin, const double cellSizeFactor = 1.0, const size_t sortingThreshold = 8,
                        LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell)
      : CellBasedParticleContainer<ParticleCellType>(boxMin, boxMax, cutoff, skin, sortingThreshold),
        _cellBlock(this->_cells, boxMin, boxMax, cutoff + skin, cellSizeFactor),
        _loadEstimator(loadEstimator) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::linkedCellsReferences; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    _cellBlock.reserve(numParticles + numParticlesHaloEstimate);
  }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const Particle_T &p) override {
    addParticleLock.lock();
    _particleList.push_back(p);
    updateDirtyParticleReferences();
    addParticleLock.unlock();
  }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const Particle_T &haloParticle) override {
    addParticleLock.lock();
    _particleList.push_back(haloParticle);
    updateDirtyParticleReferences();
    addParticleLock.unlock();
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle_T &haloParticle) override {
    auto cells = _cellBlock.getNearbyHaloCells(haloParticle.getR(), this->getVerletSkin());
    for (auto cellptr : cells) {
      bool updated = internal::checkParticleInCellAndUpdateByID(*cellptr, haloParticle);
      if (updated) {
        return true;
      }
    }
    AutoPasLog(TRACE, "UpdateHaloParticle was not able to update particle: {}", haloParticle.toString());
    return false;
  }

  /**
   * Deletes all particles from the container.
   */
  void deleteAllParticles() override {
    _particleList.clearAllParticles();
    CellBasedParticleContainer<ParticleCellType>::deleteAllParticles();
  }

  /**
   * @copydoc ParticleContainerInterface::deleteHaloParticles()
   */
  void deleteHaloParticles() override {
    _particleList.clearHaloParticles();
    _cellBlock.clearHaloCells();
  }

  /**
   * Generates the load estimation function depending on _loadEstimator.
   * @return load estimator function object.
   */
  BalancedTraversal::EstimatorFunction getLoadEstimatorFunction() {
    // (Explicit) static cast required for Apple Clang (last tested version: 17.0.0)
    switch (static_cast<LoadEstimatorOption::Value>(this->_loadEstimator)) {
      case LoadEstimatorOption::squaredParticlesPerCell: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          return loadEstimators::squaredParticlesPerCell(this->_cells, cellsPerDimension, lowerCorner, upperCorner);
        };
      }
      case LoadEstimatorOption::none:
        [[fallthrough]];
      default: {
        return
            [&](const std::array<unsigned long, 3> &cellsPerDimension, const std::array<unsigned long, 3> &lowerCorner,
                const std::array<unsigned long, 3> &upperCorner) { return 1; };
      }
    }
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override { updateDirtyParticleReferences(); }

  /**
   * Updates all the References in the cells that are out of date.
   * @note Removes all dummy particles.
   */
  void updateDirtyParticleReferences() {
    if (_particleList.isDirty()) {
      if (_particleList.needsRebuild()) {
        for (auto &cell : this->_cells) {
          cell.clear();
        }
      }

      for (auto it = _particleList.beginDirty(); it < _particleList.endDirty(); it++) {
        if (it->isDummy()) {
          continue;
        }
        ParticleCellType &cell = _cellBlock.getContainingCell(it->getR());
        auto address = &(*it);
        cell.addParticleReference(address);
      }
      _particleList.markAsClean();
    }
  }

  void computeInteractions(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *traversalInterface = dynamic_cast<LCTraversalInterface *>(traversal);
    auto *cellPairTraversal = dynamic_cast<CellTraversal<ParticleCellType> *>(traversal);
    if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
      balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
    }
    if (traversalInterface && cellPairTraversal) {
      cellPairTraversal->setSortingThreshold(this->_sortingThreshold);
      cellPairTraversal->setCellsToTraverse(this->_cells);
    } else {
      utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in LinkedCellsReferences::computeInteractions. TraversalID: {}",
          traversal->getTraversalType());
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                             IteratorBehavior iteratorBehavior,
                                                             const std::array<double, 3> &boxMin,
                                                             const std::array<double, 3> &boxMax) const override {
    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }
  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
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
  std::tuple<const Particle_T *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
                                                                 IteratorBehavior iteratorBehavior,
                                                                 const std::array<double, 3> &boxMin,
                                                                 const std::array<double, 3> &boxMax) const {
    using namespace autopas::utils::ArrayMath::literals;

    std::array<double, 3> boxMinWithSafetyMargin = boxMin;
    std::array<double, 3> boxMaxWithSafetyMargin = boxMax;
    if constexpr (regionIter) {
      // We extend the search box for cells here since particles might have moved
      boxMinWithSafetyMargin -= this->getVerletSkin();
      boxMaxWithSafetyMargin += this->getVerletSkin();
    }

    // first and last relevant cell index
    const auto [startCellIndex, endCellIndex] = [&]() -> std::tuple<size_t, size_t> {
      if constexpr (regionIter) {
        // We extend the search box for cells here since particles might have moved
        return {_cellBlock.get1DIndexOfPosition(boxMinWithSafetyMargin),
                _cellBlock.get1DIndexOfPosition(boxMaxWithSafetyMargin)};
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
          advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax,
                                             boxMinWithSafetyMargin, boxMaxWithSafetyMargin, endCellIndex);
    }

    // shortcut if the given index doesn't exist
    if (cellIndex > endCellIndex) {
      return {nullptr, 0, 0};
    }
    const Particle_T *retPtr = &this->_cells[cellIndex][particleIndex];

    return {retPtr, cellIndex, particleIndex};
  }

  /**
   * @copydoc ParticleContainerInterface::deleteParticle(Particle_T &particle)
   */
  bool deleteParticle(Particle_T &particle) override {
    // This function doesn't actually delete anything as it would mess up the reference structure.
    internal::markParticleAsDeleted(particle);
    return false;
  }

  /**
   * @copydoc ParticleContainerInterface::deleteParticle(size_t cellIndex, size_t particleIndex)
   */
  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {
    // This function doesn't actually delete anything as it would mess up the reference structure.
    internal::markParticleAsDeleted(this->_cells[cellIndex][particleIndex]);
    return false;
  }

  std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }

    std::vector<Particle_T> invalidParticles;

    // for exception handling in parallel region
    bool exceptionCaught{false};
    std::string exceptionMsg{""};

    AUTOPAS_OPENMP(parallel) {
      // private for each thread!
      std::vector<Particle_T> myInvalidParticles, myInvalidNotOwnedParticles;
      AUTOPAS_OPENMP(for)
      for (size_t cellId = 0; cellId < this->getCells().size(); ++cellId) {
        // Delete dummy particles of each cell.
        this->getCells()[cellId].deleteDummyParticles();

        // if empty
        if (this->getCells()[cellId].isEmpty()) continue;

        auto [cellLowerCorner, cellUpperCorner] = this->getCellBlock().getCellBoundingBox(cellId);

        auto &particleVec = this->getCells()[cellId]._particles;
        for (auto pIter = particleVec.begin(); pIter != particleVec.end();) {
          if ((*pIter)->isOwned() and utils::notInBox((*pIter)->getR(), cellLowerCorner, cellUpperCorner)) {
            myInvalidParticles.push_back(**pIter);

            // multi layer swap-delete
            // we copy the particle behind the last pointer in the cell over the particle we want to delete
            **pIter = *particleVec.back();
            // since we will not actually delete the particle we just copied mark it as dummy.
            particleVec.back()->setOwnershipState(OwnershipState::dummy);
            // then we pop the last pointer from the cell.
            particleVec.pop_back();

          } else {
            ++pIter;
          }
        }
      }
      // implicit barrier here
      // the barrier is needed because iterators are not thread safe w.r.t. addParticle()

      // this loop is executed for every thread and thus parallel. Don't use #pragma omp for here!
      // addParticle might throw. Set exceptionCaught to true, so we can handle that after the OpenMP region
      try {
        for (auto &&p : myInvalidParticles) {
          // if not in halo
          if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
            this->template addParticle<false>(p);
          } else {
            myInvalidNotOwnedParticles.push_back(p);
          }
        }
      } catch (const std::exception &e) {
        exceptionCaught = true;
        AUTOPAS_OPENMP(critical)
        exceptionMsg.append(e.what());
      }
      AUTOPAS_OPENMP(critical) {
        // merge private vectors to global one.
        invalidParticles.insert(invalidParticles.end(), myInvalidNotOwnedParticles.begin(),
                                myInvalidNotOwnedParticles.end());
      }
    }

    if (exceptionCaught) {
      autopas::utils::ExceptionHandler::exception(exceptionMsg);
    }

    // we have to remove halo particles after the above for-loop since removing halo particles changes the underlying
    // datastructure (_particleList) which invalidates the references in the cells.
    // Note: removing halo particles before the loop and do a call to updateDirtyParticleReferences() right after won't
    // work since there are owned particles in _particleList which fall into the halo area (they are not yet filtered
    // out by the above loop) and can therefore not added to halo cells during particle insert in
    // updateDirtyParticleReferences()
    this->deleteHaloParticles();

    _particleList.deleteDummyParticles();
    updateDirtyParticleReferences();
    return invalidParticles;
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    return TraversalSelectorInfo(this->getCellBlock().getCellsPerDimensionWithHalo(), this->getInteractionLength(),
                                 this->getCellBlock().getCellLength(), 0);
  }

  /**
   * @copydoc ParticleContainerInterface::begin()
   */
  ContainerIterator<Particle_T, true, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc ParticleContainerInterface::begin()
   */
  ContainerIterator<Particle_T, false, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors =
          nullptr) const override {
    return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHaloOrDummy) {
    if (behavior == IteratorBehavior::ownedOrHaloOrDummy) {
      // iterate over all particles, so execute directly on particle vector
      _particleList.forEach(forEachLambda);
    } else {
      for (size_t index = 0; index < getCells().size(); index++) {
        if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
          getCells()[index].forEach(forEachLambda, behavior);
        }
      }
    }
  }

  /**
   * @copydoc LinkedCells::reduce()
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHaloOrDummy) {
    if (behavior == IteratorBehavior::ownedOrHaloOrDummy) {
      // iterate over all particles, so execute directly on particle vector
      _particleList.reduce(reduceLambda, result);
    } else {
      for (size_t index = 0; index < getCells().size(); index++) {
        if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
          getCells()[index].reduce(reduceLambda, result, behavior);
        }
      }
    }
  }

  /**
   * @copydoc ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle_T, true, true>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle_T, false, true>::ParticleVecType *additionalVectors =
          nullptr) const override {
    return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc LinkedCells::forEachInRegion()
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    using namespace autopas::utils::ArrayMath::literals;

    // We increase the search region by skin, as particles can move over cell borders.
    const auto startIndex3D = this->_cellBlock.get3DIndexOfPosition(lowerCorner - this->getVerletSkin());
    const auto stopIndex3D = this->_cellBlock.get3DIndexOfPosition(higherCorner + this->getVerletSkin());

    size_t numCellsOfInterest = (stopIndex3D[0] - startIndex3D[0] + 1) * (stopIndex3D[1] - startIndex3D[1] + 1) *
                                (stopIndex3D[2] - startIndex3D[2] + 1);
    std::vector<size_t> cellsOfInterest(numCellsOfInterest);

    int i = 0;
    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          cellsOfInterest[i++] =
              utils::ThreeDimensionalMapping::threeToOneD({x, y, z}, this->_cellBlock.getCellsPerDimensionWithHalo());
        }
      }
    }

    for (size_t index : cellsOfInterest) {
      if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
        getCells()[index].forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      }
    }
  }

  /**
   * @copydoc LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    using namespace autopas::utils::ArrayMath::literals;

    // We increase the search region by skin, as particles can move over cell borders.
    const auto startIndex3D = this->_cellBlock.get3DIndexOfPosition(lowerCorner - this->getVerletSkin());
    const auto stopIndex3D = this->_cellBlock.get3DIndexOfPosition(higherCorner + this->getVerletSkin());

    size_t numCellsOfInterest = (stopIndex3D[0] - startIndex3D[0] + 1) * (stopIndex3D[1] - startIndex3D[1] + 1) *
                                (stopIndex3D[2] - startIndex3D[2] + 1);
    std::vector<size_t> cellsOfInterest(numCellsOfInterest);

    int i = 0;
    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          cellsOfInterest[i++] =
              utils::ThreeDimensionalMapping::threeToOneD({x, y, z}, this->_cellBlock.getCellsPerDimensionWithHalo());
        }
      }
    }

    for (size_t index : cellsOfInterest) {
      if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
        getCells()[index].reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
      }
    }
  }

  /**
   * Get the cell block, not supposed to be used except by verlet lists
   * @return the cell block
   */
  internal::CellBlock3D<ParticleCellType> &getCellBlock() { return _cellBlock; }

  /**
   * @copydoc getCellBlock()
   * @note const version
   */
  const internal::CellBlock3D<ParticleCellType> &getCellBlock() const { return _cellBlock; }

  /**
   * Returns reference to the data of LinkedCellsReferences
   * @return the data
   */
  std::vector<ParticleCellType> &getCells() { return this->_cells; }

 protected:
  /**
   * Given a pair of cell-/particleIndex and iterator restrictions either returns the next indices that match these
   * restrictions or indices that are out of bounds (e.g. cellIndex >= cells.size())
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin The actual search box min
   * @param boxMax The actual search box max
   * @param boxMinWithSafetyMargin Search box min that includes a surrounding of skin
   * @param boxMaxWithSafetyMargin Search box max that includes a surrounding of skin
   * @param endCellIndex Last relevant cell index
   * @return tuple<cellIndex, particleIndex>
   */
  template <bool regionIter>
  std::tuple<size_t, size_t> advanceIteratorIndices(
      size_t cellIndex, size_t particleIndex, IteratorBehavior iteratorBehavior, const std::array<double, 3> &boxMin,
      const std::array<double, 3> &boxMax, const std::array<double, 3> &boxMinWithSafetyMargin,
      const std::array<double, 3> &boxMaxWithSafetyMargin, size_t endCellIndex) const {
    // Finding the indices for the next particle
    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();

    // helper function to determine if the cell can even contain particles of interest to the iterator
    auto cellIsRelevant = [&]() -> bool {
      bool isRelevant =
          // behavior matches possible particle ownership
          (iteratorBehavior & IteratorBehavior::owned and _cellBlock.cellCanContainOwnedParticles(cellIndex)) or
          (iteratorBehavior & IteratorBehavior::halo and _cellBlock.cellCanContainHaloParticles(cellIndex));
      if constexpr (regionIter) {
        // short circuit if already false
        if (isRelevant) {
          // is the cell in the region?
          const auto [cellLowCorner, cellHighCorner] = _cellBlock.getCellBoundingBox(cellIndex);
          isRelevant =
              utils::boxesOverlap(cellLowCorner, cellHighCorner, boxMinWithSafetyMargin, boxMaxWithSafetyMargin);
        }
      }
      return isRelevant;
    };

    do {
      // advance to the next particle
      ++particleIndex;
      // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.

      // If cell has wrong type, or there are no more particles in this cell jump to the next
      while (not cellIsRelevant() or particleIndex >= this->_cells[cellIndex].size()) {
        // TODO: can this jump be done more efficient if behavior is only halo or owned?
        // TODO: can this jump be done more efficient for region iters if the cell is outside the region?
        cellIndex += stride;
        particleIndex = 0;

        // If we notice that there is nothing else to look at set invalid values, so we get a nullptr next time and
        // break.
        if (cellIndex > endCellIndex) {
          return {std::numeric_limits<decltype(cellIndex)>::max(), std::numeric_limits<decltype(particleIndex)>::max()};
        }
      }
    } while (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
        this->_cells[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax));

    // the indices returned at this point should always be valid
    return {cellIndex, particleIndex};
  }

  /**
   * object that stores the actual Particles and keeps track of the references.
   */
  ParticleVector<Particle_T> _particleList;
  /**
   * object to manage the block of cells.
   */
  internal::CellBlock3D<ParticleCellType> _cellBlock;
  /**
   * load estimation algorithm for balanced traversals.
   */
  LoadEstimatorOption _loadEstimator;
  /**
   * Workaround for adding particles in parallel -> https://github.com/AutoPas/AutoPas/issues/555
   */
  AutoPasLock addParticleLock;
};

}  // namespace autopas
