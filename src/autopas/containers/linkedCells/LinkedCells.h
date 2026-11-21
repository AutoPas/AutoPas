/**
 * @file LinkedCells.h
 *
 * @author tchipevn
 * @date 17.02.2018
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/containers/cellTraversals/BalancedTraversal.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCTraversalInterface.h"
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
class LinkedCells : public CellBasedParticleContainer<FullParticleCell<Particle_T>> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCellType = FullParticleCell<Particle_T>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle_T;

  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param aosSortingThresholdFallback number of particles in two cells from which AoS sorting should be performed.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
   * @param soaSortingThresholdFallback Sum of the SoA buffer sizes of two cells from which SoA sorting should be
   * enabled. Default comes from the LJFunctorHWY Benchmarks.
   * @param loadEstimator the load estimation algorithm for balanced traversals.
   * By default all applicable traversals are allowed.
   */
  LinkedCells(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
              const double skin, const double cellSizeFactor = 1.0, const size_t aosSortingThresholdFallback = 8,
              const size_t soaSortingThresholdFallback = 100,
              LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell)
      : CellBasedParticleContainer<ParticleCellType>(boxMin, boxMax, cutoff, skin, aosSortingThresholdFallback,
                                                     soaSortingThresholdFallback),
        _cellBlock(this->_cells, boxMin, boxMax, cutoff + skin, cellSizeFactor),
        _loadEstimator(loadEstimator) {}

  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::linkedCells; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    _cellBlock.reserve(numParticles + numParticlesHaloEstimate);
  }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const Particle_T &p) override {
    ParticleCellType &cell = _cellBlock.getContainingCell(p.getR());
    cell.addParticle(p);
  }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const Particle_T &haloParticle) override {
    ParticleCellType &cell = _cellBlock.getContainingCell(haloParticle.getR());
    cell.addParticle(haloParticle);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle_T &haloParticle) override {
    const auto skinHalf = 0.5 * this->getVerletSkin();
    const bool updated = _cellBlock.forEachNearbyHaloCell(haloParticle.getR(), skinHalf, [&](auto &cell) {
      return internal::checkParticleInCellAndUpdateByIDAndPosition(cell, haloParticle, skinHalf);
    });
    if (updated) {
      return true;
    }
    AutoPasLog(TRACE, "UpdateHaloParticle was not able to update particle: {}", haloParticle.toString());
    return false;
  }

  /**
   * @copydoc ParticleContainerInterface::deleteHaloParticles()
   */
  void deleteHaloParticles() override { _cellBlock.clearHaloCells(); }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
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

  void computeInteractions(TraversalInterface *traversal) override {
    prepareTraversal(traversal);

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }

    this->deleteHaloParticles();

    std::vector<Particle_T> invalidParticles;
    AUTOPAS_OPENMP(parallel) {
      // private for each thread!
      std::vector<Particle_T> myInvalidParticles{}, myInvalidNotOwnedParticles{};
      // TODO: needs smarter heuristic than this.
      myInvalidParticles.reserve(128);
      myInvalidNotOwnedParticles.reserve(128);
      AUTOPAS_OPENMP(for)
      for (size_t cellId = 0; cellId < this->getCells().size(); ++cellId) {
        // Delete dummy particles of each cell.
        this->getCells()[cellId].deleteDummyParticles();

        // if empty
        if (this->getCells()[cellId].isEmpty()) continue;

        const auto [cellLowerCorner, cellUpperCorner] = this->getCellBlock().getCellBoundingBox(cellId);

        auto &particleVec = this->getCells()[cellId]._particles;
        for (auto pIter = particleVec.begin(); pIter < particleVec.end();) {
          // if not in cell
          if (utils::notInBox(pIter->getR(), cellLowerCorner, cellUpperCorner)) {
            myInvalidParticles.push_back(*pIter);
            // swap-delete
            *pIter = particleVec.back();
            particleVec.pop_back();
          } else {
            ++pIter;
          }
        }
      }
      // implicit barrier here
      // the barrier is needed because iterators are not threadsafe w.r.t. addParticle()

      // this loop is executed for every thread and thus parallel. Don't use #pragma omp for here!
      for (auto &&p : myInvalidParticles) {
        // if not in halo
        if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
          this->template addParticle<false>(p);
        } else {
          myInvalidNotOwnedParticles.push_back(p);
        }
      }
      AUTOPAS_OPENMP(critical) {
        // merge private vectors to global one.
        invalidParticles.insert(invalidParticles.end(), myInvalidNotOwnedParticles.begin(),
                                myInvalidNotOwnedParticles.end());
      }
    }
    return invalidParticles;
  }

  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    return TraversalSelectorInfo(this->getCellBlock().getCellsPerDimensionWithHalo(), this->getInteractionLength(),
                                 this->getCellBlock().getCellLength(), 0);
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

    if (cellIndex >= this->_cells.size()) {
      return {nullptr, 0, 0};
    }

    const bool isStart = (cellIndex == 0 and particleIndex == 0);

    // Fast path: when advancing within an already-validated cell, scan its remaining particles directly without
    // recomputing cell bounds or 3D indices.
    if (not isStart) {
      const auto &cell = this->_cells[cellIndex];
      for (; particleIndex < cell.size(); ++particleIndex) {
        if (containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
                cell[particleIndex], iteratorBehavior, boxMin, boxMax)) {
          return {&cell[particleIndex], cellIndex, particleIndex};
        }
      }
    }

    const bool includeOwned = static_cast<bool>(iteratorBehavior & IteratorBehavior::owned);
    const bool includeHalo = static_cast<bool>(iteratorBehavior & IteratorBehavior::halo);
    if (not includeOwned and not includeHalo) {
      return {nullptr, 0, 0};
    }

    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();
    const size_t threadOffset = (iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num();

    if constexpr (regionIter) {
      // Extend the search box for cells by verletSkin since particles might have drifted since the last rebuild.
      const auto boxMinWithSafetyMargin = boxMin - this->getVerletSkin();
      const auto boxMaxWithSafetyMargin = boxMax + this->getVerletSkin();

      auto minIdx3D = _cellBlock.get3DIndexOfPosition(boxMinWithSafetyMargin);
      auto maxIdx3D = _cellBlock.get3DIndexOfPosition(boxMaxWithSafetyMargin);

      const auto &dims = _cellBlock.getCellsPerDimensionWithHalo();
      const auto cpl = _cellBlock.getCellsPerInteractionLength();

      // If only owned particles are requested, clamp the 3D sub-box directly to owned cells.
      if (not includeHalo) {
        for (size_t d = 0; d < 3; ++d) {
          minIdx3D[d] = std::max(minIdx3D[d], cpl);
          maxIdx3D[d] = std::min(maxIdx3D[d], dims[d] - cpl - 1);
          if (minIdx3D[d] > maxIdx3D[d]) {
            return {nullptr, 0, 0};
          }
        }
      }

      const std::array<size_t, 3> subDims{
          maxIdx3D[0] - minIdx3D[0] + 1,
          maxIdx3D[1] - minIdx3D[1] + 1,
          maxIdx3D[2] - minIdx3D[2] + 1,
      };
      const size_t totalSubCells = subDims[0] * subDims[1] * subDims[2];

      size_t subIdx{};
      if (isStart) {
        subIdx = threadOffset;
      } else {
        const auto currentIdx3D = utils::ThreeDimensionalMapping::oneToThreeD(cellIndex, dims);
        subIdx = ((currentIdx3D[2] - minIdx3D[2]) * subDims[1] + (currentIdx3D[1] - minIdx3D[1])) * subDims[0] +
                 (currentIdx3D[0] - minIdx3D[0]) + stride;
      }

      const bool haloOnly = includeHalo and not includeOwned;
      for (; subIdx < totalSubCells; subIdx += stride) {
        const auto local3D = utils::ThreeDimensionalMapping::oneToThreeD(subIdx, subDims);
        const std::array<size_t, 3> cell3D{
            minIdx3D[0] + local3D[0],
            minIdx3D[1] + local3D[1],
            minIdx3D[2] + local3D[2],
        };
        if (haloOnly) {
          const bool isHaloCell = cell3D[0] < cpl or cell3D[0] >= dims[0] - cpl or cell3D[1] < cpl or
                                  cell3D[1] >= dims[1] - cpl or cell3D[2] < cpl or cell3D[2] >= dims[2] - cpl;
          if (not isHaloCell) {
            continue;
          }
        }
        const size_t nextCellIndex = utils::ThreeDimensionalMapping::threeToOneD(cell3D, dims);
        const auto &cell = this->_cells[nextCellIndex];
        for (size_t pIdx = 0; pIdx < cell.size(); ++pIdx) {
          if (containerIteratorUtils::particleFulfillsIteratorRequirements<true>(cell[pIdx], iteratorBehavior, boxMin,
                                                                                 boxMax)) {
            return {&cell[pIdx], nextCellIndex, pIdx};
          }
        }
      }
      return {nullptr, 0, 0};
    } else {
      const auto [startCellIndex, endCellIndex] =
          includeHalo
              ? std::tuple<size_t, size_t>{0, this->_cells.size() - 1}
              : std::tuple<size_t, size_t>{_cellBlock.getFirstOwnedCellIndex(), _cellBlock.getLastOwnedCellIndex()};

      cellIndex = isStart ? (startCellIndex + threadOffset) : (cellIndex + stride);
      const bool checkCellOwnership = not(includeOwned and includeHalo);
      for (; cellIndex <= endCellIndex; cellIndex += stride) {
        if (checkCellOwnership) {
          const bool isRelevant = (includeOwned and _cellBlock.cellCanContainOwnedParticles(cellIndex)) or
                                  (includeHalo and _cellBlock.cellCanContainHaloParticles(cellIndex));
          if (not isRelevant) {
            continue;
          }
        }
        const auto &cell = this->_cells[cellIndex];
        for (size_t pIdx = 0; pIdx < cell.size(); ++pIdx) {
          if (containerIteratorUtils::particleFulfillsIteratorRequirements<false>(cell[pIdx], iteratorBehavior, boxMin,
                                                                                  boxMax)) {
            return {&cell[pIdx], cellIndex, pIdx};
          }
        }
      }
      return {nullptr, 0, 0};
    }
  }

  /**
   * @copydoc ParticleContainerInterface::deleteParticle(Particle_T &particle)
   */
  void deleteParticle(Particle_T &particle) override { internal::markParticleAsDeleted(particle); }

  /**
   * @copydoc ParticleContainerInterface::deleteParticle(size_t cellIndex, size_t particleIndex)
   */
  void deleteParticle(size_t cellIndex, size_t particleIndex) override {
    // This function doesn't actually delete anything as it would mess up the references in neighbor lists that are
    // built on top of this container.
    internal::markParticleAsDeleted(this->_cells[cellIndex][particleIndex]);
  }

  /**
   * @copydoc ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      utils::optRef<typename ContainerIterator<Particle_T, true, false>::ParticleVecType> additionalVectors =
          std::nullopt) override {
    return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      utils::optRef<typename ContainerIterator<Particle_T, false, false>::ParticleVecType> additionalVectors =
          std::nullopt) const override {
    return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
  }

  /**
   * Execute code on all particles in this container as defined by a lambda function.
   * @tparam Lambda (ParticleType &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    if (behavior == IteratorBehavior::ownedOrHaloOrDummy) {
      for (auto &cell : getCells()) {
        cell.forEach(forEachLambda);
      }
    } else {
      for (size_t index = 0; index < getCells().size(); index++) {
        if (not _cellBlock.ignoreCellForIteration(index, behavior)) {
          getCells()[index].forEach(forEachLambda, behavior);
        }
      }
    }
  }

  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (ParticleType p, A initialValue) -> void
   * @tparam A type of particle attribute to be reduced
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownedOrHalo
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    if (behavior == IteratorBehavior::ownedOrHaloOrDummy) {
      for (auto &cell : getCells()) {
        cell.reduce(reduceLambda, result);
      }
    } else {
      for (size_t index = 0; index < getCells().size(); index++) {
        if (not _cellBlock.ignoreCellForIteration(index, behavior)) {
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
      utils::optRef<typename ContainerIterator<Particle_T, true, true>::ParticleVecType> additionalVectors =
          std::nullopt) override {
    return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      utils::optRef<typename ContainerIterator<Particle_T, false, true>::ParticleVecType> additionalVectors =
          std::nullopt) const override {
    return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * Execute code on all particles in this container in a certain region as defined by a lambda function.
   * @tparam Lambda (ParticleType &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    using namespace autopas::utils::ArrayMath::literals;

    const auto startIndex3D = this->_cellBlock.get3DIndexOfPosition(lowerCorner - this->getVerletSkin());
    const auto stopIndex3D = this->_cellBlock.get3DIndexOfPosition(higherCorner + this->getVerletSkin());

    const size_t numCellsOfInterest = (stopIndex3D[0] - startIndex3D[0] + 1) * (stopIndex3D[1] - startIndex3D[1] + 1) *
                                      (stopIndex3D[2] - startIndex3D[2] + 1);
    std::vector<size_t> cellsOfInterest;
    cellsOfInterest.reserve(numCellsOfInterest);

    const auto &cellsPerDimensionWithHalo = this->_cellBlock.getCellsPerDimensionWithHalo();

    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          cellsOfInterest.push_back(utils::ThreeDimensionalMapping::threeToOneD({x, y, z}, cellsPerDimensionWithHalo));
        }
      }
    }

    for (auto cellIndex : cellsOfInterest) {
      if (not _cellBlock.ignoreCellForIteration(cellIndex, behavior)) {
        getCells()[cellIndex].forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      }
    }
  }

  /**
   * Execute code on all particles in this container in a certain region as defined by a lambda function.
   * @tparam Lambda (ParticleType &p, A &result) -> void
   * @tparam A type of reduction Value
   * @param reduceLambda code to be executed on all particles
   * @param result reference to starting and final value for reduction
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    using namespace autopas::utils::ArrayMath::literals;
    const auto startIndex3D = this->_cellBlock.get3DIndexOfPosition(lowerCorner - this->getVerletSkin());
    const auto stopIndex3D = this->_cellBlock.get3DIndexOfPosition(higherCorner + this->getVerletSkin());

    const size_t numCellsOfInterest = (stopIndex3D[0] - startIndex3D[0] + 1) * (stopIndex3D[1] - startIndex3D[1] + 1) *
                                      (stopIndex3D[2] - startIndex3D[2] + 1);
    std::vector<size_t> cellsOfInterest;
    cellsOfInterest.reserve(numCellsOfInterest);

    const auto &cellsPerDimensionWithHalo = this->_cellBlock.getCellsPerDimensionWithHalo();

    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          cellsOfInterest.push_back(utils::ThreeDimensionalMapping::threeToOneD({x, y, z}, cellsPerDimensionWithHalo));
        }
      }
    }

    for (auto cellIndex : cellsOfInterest) {
      if (not _cellBlock.ignoreCellForIteration(cellIndex, behavior)) {
        getCells()[cellIndex].reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
      }
    }
  }

  /**
   * Get the cell block, not supposed to be used except by verlet lists
   * @return the cell block
   */
  internal::CellBlock3D<ParticleCellType> &getCellBlock() { return _cellBlock; }

  /**
   * @copydoc autopas::LinkedCells::getCellBlock()
   * @note const version
   */
  const internal::CellBlock3D<ParticleCellType> &getCellBlock() const { return _cellBlock; }

  /**
   * Returns a non-const reference to the cell data structure.
   * @return Non-const reference.
   */
  std::vector<ParticleCellType> &getCells() { return this->_cells; }

 protected:
  /**
   * Checks if a given traversal is allowed for LinkedCells and sets it up for the force interactions.
   * @tparam Traversal Traversal type. E.g. pairwise, triwise
   * @param traversal
   */
  template <typename Traversal>
  void prepareTraversal(Traversal &traversal) {
    auto *traversalInterface = dynamic_cast<LCTraversalInterface *>(traversal);
    auto *cellTraversal = dynamic_cast<CellTraversal<ParticleCellType> *>(traversal);
    if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
      balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
    }
    if (traversalInterface && cellTraversal) {
      cellTraversal->setAoSSortingThresholds(*this->_aosSortingThresholds);
      cellTraversal->setSoASortingThresholds(*this->_soaSortingThresholds);
      cellTraversal->setCellsToTraverse(this->_cells);
    } else {
      utils::ExceptionHandler::exception(
          "The selected traversal is not compatible with the LinkedCells container. TraversalID: {}",
          traversal->getTraversalType());
    }
  }

  /**
   * object to manage the block of cells.
   */
  internal::CellBlock3D<ParticleCellType> _cellBlock;

  /**
   * load estimation algorithm for balanced traversals.
   */
  LoadEstimatorOption _loadEstimator;
};

}  // namespace autopas
