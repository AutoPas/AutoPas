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
#include "autopas/containers/cellTraversals/CellTraversal.h"
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
 * This class stores all owned particles in a single cell.
 * Interactions are calculated directly, such that each particle interacts with
 * every other particle.
 * Use this class only if you have a very small amount of particles at hand.
 *
 * The surrounding halo volume is decomposed into six halo cells. This is done, so that cell sorting can be used for
 * interactions with halo particles. The volumes are not equal in size, which should not matter since we only use a
 * sequential traversal (so far). The image below shows how the halo volume is split up:
 *
 * \image html DSHalos.svg "Halo Volume Segmentation" width=75%
 *
 * The image above depicts the segmentation of the surrounding halo volume into six cells, one cell on every side,
 * where the halo cell centers always align with the owned cell center.
 *
 * @tparam Particle_T Particle type that is used with this container.
 */
template <class Particle_T>
class DirectSum : public CellBasedParticleContainer<FullParticleCell<Particle_T>> {
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
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param sortingThreshold
   */
  DirectSum(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double cutoff, double skin,
            const size_t sortingThreshold)
      : CellBasedParticleContainer<ParticleCellType>(boxMin, boxMax, cutoff, skin, sortingThreshold),
        _cellBorderFlagManager() {
    using namespace autopas::utils::ArrayMath::literals;
    // 1 owned and 6 halo cells
    this->_cells.resize(7);
    this->_cells[0].setPossibleParticleOwnerships(OwnershipState::owned);
    std::for_each(++this->_cells.begin(), this->_cells.end(),
                  [&](auto &cell) { cell.setPossibleParticleOwnerships(OwnershipState::halo); });
    auto boxLength = boxMax - boxMin;
    this->_cells[0].setCellLength(boxLength);
  }

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::directSum; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    this->getOwnedCell().reserve(numParticles);
    for (auto cellIt = ++this->_cells.begin(); cellIt != this->_cells.end(); cellIt++) {
      cellIt->reserve(numParticlesHaloEstimate);
    }
  };

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const Particle_T &p) override { getOwnedCell().addParticle(p); }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const Particle_T &haloParticle) override {
    const auto boxMax = this->getBoxMax();
    const auto boxMin = this->getBoxMin();
    const auto pos = haloParticle.getR();

    for (size_t dim = 0; dim < 3; ++dim) {
      if (pos[dim] < boxMin[dim]) {
        this->_cells[2 * dim + 1].addParticle(haloParticle);
        return;
      } else if (pos[dim] >= boxMax[dim]) {
        this->_cells[2 * dim + 2].addParticle(haloParticle);
        return;
      }
    }
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle_T &haloParticle) override {
    const auto boxMax = this->getBoxMax();
    const auto boxMin = this->getBoxMin();
    const auto pos = haloParticle.getR();
    const auto skinHalf = 0.5 * this->getVerletSkin();

    // Look for the particle in halo cells that are within half the skin distance of its position
    for (size_t dim = 0; dim < 3; ++dim) {
      if (pos[dim] < boxMin[dim] + skinHalf) {
        if (internal::checkParticleInCellAndUpdateByIDAndPosition(this->_cells[2 * dim + 1], haloParticle, skinHalf)) {
          return true;
        }
      } else if (pos[dim] >= boxMax[dim] - skinHalf) {
        if (internal::checkParticleInCellAndUpdateByIDAndPosition(this->_cells[2 * dim + 2], haloParticle, skinHalf)) {
          return true;
        }
      }
    }
    return false;
  }

  void deleteHaloParticles() override {
    for (auto cellIt = ++this->_cells.begin(); cellIt != this->_cells.end(); cellIt++) {
      cellIt->clear();
    }
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
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
    // first we delete halo particles, as we don't want them here.
    deleteHaloParticles();
    getOwnedCell().deleteDummyParticles();

    std::vector<Particle_T> invalidParticles{};
    auto &particleVec = getOwnedCell()._particles;
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
    using namespace autopas::utils::ArrayMath::literals;

    // direct sum consists of seven cells (owned + two halo cells in each dimension)
    return TraversalSelectorInfo(
        // DS container can be viewed as a 3x3x3 grid with some halo cells being combined.
        {3, 3, 3},
        // intentionally use cutoff here, as ds_sequential should be using the cutoff.
        this->getCutoff(), this->getBoxMax() - this->getBoxMin(), 0);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors =
          nullptr) const override {
    return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    if (behavior & IteratorBehavior::owned) {
      getOwnedCell().forEach(forEachLambda);
    }
    if (behavior & IteratorBehavior::halo) {
      for (auto cellIt = ++this->_cells.begin(); cellIt != this->_cells.end(); cellIt++) {
        cellIt->forEach(forEachLambda);
      }
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
      getOwnedCell().reduce(reduceLambda, result);
    }
    if (behavior & IteratorBehavior::halo) {
      for (auto cellIt = ++this->_cells.begin(); cellIt != this->_cells.end(); cellIt++) {
        cellIt->reduce(reduceLambda, result);
      }
    }
    // sanity check
    if (not(behavior & IteratorBehavior::ownedOrHalo)) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle_T, true, true>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
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
    if (behavior & IteratorBehavior::owned) {
      getOwnedCell().forEach(forEachLambda, lowerCorner, higherCorner, behavior);
    }
    if (behavior & IteratorBehavior::halo) {
      for (auto cellIt = ++this->_cells.begin(); cellIt != this->_cells.end(); cellIt++) {
        cellIt->forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      }
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
      getOwnedCell().reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
    }
    if (behavior & IteratorBehavior::halo) {
      for (auto cellIt = ++this->_cells.begin(); cellIt != this->_cells.end(); cellIt++) {
        cellIt->reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
      }
    }
    // sanity check
    if (not(behavior & IteratorBehavior::ownedOrHalo)) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
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
   * @copydoc autopas::ParticleContainerInterface::deleteParticle()
   */
  bool deleteParticle(Particle_T &particle) override {
    // swap-delete helper function
    auto swapDelFromCell = [&](auto &particleCell) -> bool {
      auto &particleVec = particleCell._particles;
      const bool isRearParticle = &particle == &particleVec.back();
      particle = particleVec.back();
      particleVec.pop_back();
      return isRearParticle;
    };

    // deduce into which vector the reference points
    if (particle.isOwned()) {
      return swapDelFromCell(getOwnedCell());
    } else if (particle.isHalo()) {
      const auto boxMin = this->getBoxMin();
      const auto boxMax = this->getBoxMax();
      const auto skinHalf = 0.5 * this->getVerletSkin();
      const auto pos = particle.getR();

      // Look for the particle in halo cells that are within half the skinHalf distance of its position
      for (size_t dim = 0; dim < 3; ++dim) {
        if (pos[dim] < boxMin[dim] + skinHalf) {
          if (swapDelFromCell(this->_cells[2 * dim + 1])) {
            return true;
          }
        } else if (pos[dim] >= boxMax[dim] - skinHalf) {
          if (swapDelFromCell(this->_cells[2 * dim + 2])) {
            return true;
          }
        }
      }
    }
    return false;
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
    [[nodiscard]] bool cellCanContainHaloParticles(index_t index1d) const override {
      return index1d >= 1 and index1d <= 6;
    }

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
  std::tuple<const Particle_T *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
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
        return {1, 6};
      }
      if constexpr (regionIter) {
        // if the region lies fully within the container only look at the owned cell
        if (utils::ArrayMath::less(this->getBoxMin(), boxMin) and utils::ArrayMath::less(boxMax, this->getBoxMax())) {
          return {0, 0};
        }
      }
      // all cells
      return {0, 6};
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
    if (particleIndex >= this->_cells[cellIndex].size() or
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
    const Particle_T *retPtr = &this->_cells[cellIndex][particleIndex];

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
    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_preferred_num_threads();

    do {
      // advance to the next particle
      ++particleIndex;
      // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.
      while (particleIndex >= this->_cells[cellIndex].size()) {
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

  /**
   * Checks if a given traversal is allowed for DirectSum and sets it up for the force interactions.
   * @tparam Traversal Traversal type. E.g. pairwise, triwise
   * @param traversal
   */
  template <typename Traversal>
  void prepareTraversal(Traversal &traversal) {
    auto *dsTraversal = dynamic_cast<DSTraversalInterface *>(traversal);
    auto *cellTraversal = dynamic_cast<CellTraversal<ParticleCellType> *>(traversal);
    if (dsTraversal && cellTraversal) {
      cellTraversal->setSortingThreshold(this->_sortingThreshold);
      cellTraversal->setCellsToTraverse(this->_cells);
    } else {
      utils::ExceptionHandler::exception(
          "The selected traversal is not compatible with the DirectSum container. TraversalID: {}",
          traversal->getTraversalType());
    }
  }

  ParticleCellType &getOwnedCell() { return this->_cells[0]; };

  ParticleCellType &getHaloCell() { return this->_cells[1]; };
};

}  // namespace autopas
