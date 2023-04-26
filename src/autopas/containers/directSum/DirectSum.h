/**
 * @file DirectSum.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/directSum/traversals/DSTraversalInterface.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AutoPasMacros.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ParticleCellHelpers.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * This class stores particles in a single cell.
 * Interactions are calculated directly, such that each particle interacts with
 * every other particle.
 * Use this class only if you have a very small amount of particles at hand.
 * @tparam ParticleCell type of the cell that stores the particle
 */
template <class Particle>
class DirectSum : public CellBasedParticleContainer<FullParticleCell<Particle>> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = typename CellBasedParticleContainer<ParticleCell>::ParticleType;

  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skinPerTimestep
   * @param verletRebuildFrequency
   */
  DirectSum(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
            double skinPerTimestep, unsigned int verletRebuildFrequency)
      : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skinPerTimestep * verletRebuildFrequency),
        _cellBorderFlagManager() {
    this->_cells.resize(2);
  }

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::directSum; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    this->getCell().reserve(numParticles);
    this->getHaloCell().reserve(numParticlesHaloEstimate);
  };

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const ParticleType &p) override { getCell().addParticle(p); }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    ParticleType p_copy = haloParticle;
    p_copy.setOwnershipState(OwnershipState::halo);
    getHaloCell().addParticle(p_copy);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwnershipState(OwnershipState::halo);
    return internal::checkParticleInCellAndUpdateByIDAndPosition(getHaloCell(), pCopy, this->getVerletSkin());
  }

  void deleteHaloParticles() override { getHaloCell().clear(); }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
  }

  bool neighborListsAreValid() override { return true; }

  CellType getParticleCellTypeEnum() override { return CellType::FullParticleCell; }

  void iteratePairwise(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *traversalInterface = dynamic_cast<DSTraversalInterface<ParticleCell> *>(traversal);
    auto *cellPairTraversal = dynamic_cast<CellPairTraversal<ParticleCell> *>(traversal);
    if (traversalInterface && cellPairTraversal) {
      cellPairTraversal->setCellsToTraverse(this->_cells);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in DirectSum::iteratePairwise");
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return autopas::LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }
    // first we delete halo particles, as we don't want them here.
    deleteHaloParticles();
    getCell().deleteDummyParticles();

    std::vector<ParticleType> invalidParticles{};
    auto &particleVec = getCell()._particles;
    for (auto iter = particleVec.begin(); iter != particleVec.end();) {
      if (utils::notInBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        invalidParticles.push_back(*iter);
        // swap-delete
        *iter = particleVec.back();
        particleVec.pop_back();
      } else {
        ++iter;
      }
    }
    return invalidParticles;
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // direct sum technically consists of two cells (owned + halo)
    return TraversalSelectorInfo(
        {2, 0, 0},
        this->getCutoff() /*intentionally use cutoff here, as the directsumtraversal should be using the cutoff.*/,
        utils::ArrayMath::sub(this->getBoxMax(), this->getBoxMin()), 0);
  }

  [[nodiscard]] ContainerIterator<ParticleType, true, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, true, false>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<ParticleType, true, false>(*this, behavior, additionalVectors);
  }

  [[nodiscard]] ContainerIterator<ParticleType, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, false, false>::ParticleVecType *additionalVectors =
          nullptr) const override {
    return ContainerIterator<ParticleType, false, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    if (behavior & IteratorBehavior::owned) {
      getCell().forEach(forEachLambda);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().forEach(forEachLambda);
    }
    // sanity check
    if (not(behavior & IteratorBehavior::ownedOrHalo)) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  /**
   * @copydoc LinkedCells::reduce()
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior) {
    if (behavior & IteratorBehavior::owned) {
      getCell().reduce(reduceLambda, result);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().reduce(reduceLambda, result);
    }
    // sanity check
    if (not(behavior & IteratorBehavior::ownedOrHalo)) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  [[nodiscard]] ContainerIterator<ParticleType, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, true, true>::ParticleVecType *additionalVectors) override {
    return ContainerIterator<ParticleType, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  [[nodiscard]] ContainerIterator<ParticleType, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, false, true>::ParticleVecType *additionalVectors) const override {
    return ContainerIterator<ParticleType, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc LinkedCells::forEachInRegion()
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    if (behavior & IteratorBehavior::owned) {
      getCell().forEach(forEachLambda, lowerCorner, higherCorner, behavior);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().forEach(forEachLambda, lowerCorner, higherCorner, behavior);
    }
    // sanity check
    if (not(behavior & IteratorBehavior::ownedOrHalo)) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  /**
   * @copydoc LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    if (behavior & IteratorBehavior::owned) {
      getCell().reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
    }
    // sanity check
    if (not(behavior & IteratorBehavior::ownedOrHalo)) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                           IteratorBehavior iteratorBehavior,
                                                           const std::array<double, 3> &boxMin,
                                                           const std::array<double, 3> &boxMax) const override {
    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                           IteratorBehavior iteratorBehavior) const override {
    // this is not a region iter hence we stretch the bounding box to the numeric max
    constexpr std::array<double, 3> boxMin{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                                           std::numeric_limits<double>::lowest()};

    constexpr std::array<double, 3> boxMax{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                           std::numeric_limits<double>::max()};
    return getParticleImpl<false>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }

  bool deleteParticle(Particle &particle) override {
    // deduce into which vector the reference points
    auto &particleVec = particle.isOwned() ? getCell()._particles : getHaloCell()._particles;
    const bool isRearParticle = &particle == &particleVec.back();
    // swap-delete
    particle = particleVec.back();
    particleVec.pop_back();
    return not isRearParticle;
  }

  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {
    auto &particleVec = this->_cells[cellIndex]._particles;
    auto &particle = particleVec[particleIndex];
    // swap-delete
    particle = particleVec.back();
    particleVec.pop_back();
    return particleIndex < particleVec.size();
  }

 private:
  class DirectSumCellBorderAndFlagManager : public internal::CellBorderAndFlagManager {
    /**
     * the index type to access the particle cells
     */
    using index_t = std::size_t;

   public:
    [[nodiscard]] bool cellCanContainHaloParticles(index_t index1d) const override { return index1d == 1; }

    [[nodiscard]] bool cellCanContainOwnedParticles(index_t index1d) const override { return index1d == 0; }

  } _cellBorderFlagManager;

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
      // shortcuts to limit the iterator to only part of the domain.
      if (not(iteratorBehavior & IteratorBehavior::halo)) {
        // only owned region
        return {0, 0};
      }
      if (not(iteratorBehavior & IteratorBehavior::owned)) {
        // only halo region
        return {1, 1};
      }
      if constexpr (regionIter) {
        // if the region lies fully within the container only look at the owned cell
        if (utils::ArrayMath::less(this->getBoxMin(), boxMin) and utils::ArrayMath::less(boxMax, this->getBoxMax())) {
          return {0, 0};
        }
      }
      // all cells
      return {0, 1};
    }();

    // if we are at the start of an iteration determine this thread's cell index
    if (cellIndex == 0 and particleIndex == 0) {
      cellIndex =
          startCellIndex + ((iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num());
    }
    // abort if the index is out of bounds
    if (cellIndex >= this->_cells.size()) {
      return {nullptr, 0, 0};
    }
    // check the data behind the indices
    if (particleIndex >= this->_cells[cellIndex].numParticles() or
        not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
            this->_cells[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax)) {
      // either advance them to something interesting or invalidate them.
      std::tie(cellIndex, particleIndex) =
          advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
    }

    // shortcut if the given index doesn't exist
    if (cellIndex >= this->_cells.size()) {
      return {nullptr, 0, 0};
    }
    const Particle *retPtr = &this->_cells[cellIndex][particleIndex];

    return {retPtr, cellIndex, particleIndex};
  }

  /**
   * Given a pair of cell-/particleIndex and iterator restrictions either returns the next indices that match these
   * restrictions or indices that are out of bounds (e.g. cellIndex >= cells.size())
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin
   * @param boxMax
   * @return tuple<cellIndex, particleIndex>
   */
  template <bool regionIter>
  std::tuple<size_t, size_t> advanceIteratorIndices(size_t cellIndex, size_t particleIndex,
                                                    IteratorBehavior iteratorBehavior,
                                                    const std::array<double, 3> &boxMin,
                                                    const std::array<double, 3> &boxMax) const {
    // Find the indices for the next particle
    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();

    do {
      // advance to the next particle
      ++particleIndex;
      // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.
      while (particleIndex >= this->_cells[cellIndex].numParticles()) {
        cellIndex += stride;
        particleIndex = 0;
        // If there are no more reasonable cells return invalid indices.
        if (cellIndex > ((not(iteratorBehavior & IteratorBehavior::halo)) ? 0 : (this->_cells.size() - 1))) {
          return {std::numeric_limits<decltype(cellIndex)>::max(), std::numeric_limits<decltype(particleIndex)>::max()};
        }
      }
    } while (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
        this->_cells[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax));

    // the indices returned at this point should always be valid
    return {cellIndex, particleIndex};
  }

  ParticleCell &getCell() { return this->_cells[0]; };

  ParticleCell &getHaloCell() { return this->_cells[1]; };
};

}  // namespace autopas
