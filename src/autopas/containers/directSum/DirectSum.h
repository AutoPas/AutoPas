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

  [[nodiscard]] ContainerIterator<ParticleType, true> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, true>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<ParticleType, true>(*this, behavior, additionalVectors);
  }

  [[nodiscard]] ContainerIterator<ParticleType, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, false>::ParticleVecType *additionalVectors = nullptr) const override {
    return ContainerIterator<ParticleType, false>(*this, behavior, additionalVectors);
  }

  /**
   * @copydoc LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    auto forEach = [&](ParticleCell &cell) {
      for (Particle &p : cell._particles) {
        forEachLambda(p);
      }
    };

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
    auto forEach = [&](ParticleCell &cell) {
      for (Particle &p : cell._particles) {
        reduceLambda(p, result);
      }
    };

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

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                              const std::array<double, 3> &higherCorner,
                                                                              IteratorBehavior behavior) override {
    std::vector<size_t> cellsOfInterest;

    if (behavior & IteratorBehavior::owned) {
      cellsOfInterest.push_back(0);
    }
    if (behavior & IteratorBehavior::halo) {
      cellsOfInterest.push_back(1);
    }
    // sanity check
    if (cellsOfInterest.empty()) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }

    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, true>(
            &this->_cells, lowerCorner, higherCorner, std::move(cellsOfInterest), &_cellBorderFlagManager, behavior, nullptr));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) const override {
    std::vector<size_t> cellsOfInterest;

    if (behavior & IteratorBehavior::owned) {
      cellsOfInterest.push_back(0);
    }
    if (behavior & IteratorBehavior::halo) {
      cellsOfInterest.push_back(1);
    }
    // sanity check
    if (cellsOfInterest.empty()) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }

    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, false>(
            &this->_cells, lowerCorner, higherCorner, std::move(cellsOfInterest), &_cellBorderFlagManager, behavior, nullptr));
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
    // shortcut if the given index doesn't exist
    if (cellIndex >= this->_cells.size() or particleIndex >= this->_cells[cellIndex].numParticles()) {
      return {nullptr, 0, 0};
    }
    const Particle *retPtr = &this->_cells[cellIndex][particleIndex];

    // Finding the indices for the next particle
    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();
    // If we only look at halo particles and are in sequential mode, directly jump to the halo cell.
    // If there is more than one thread this will already be covered by the second directly.
    if ((iteratorBehavior & IteratorBehavior::halo) and
        ((iteratorBehavior & IteratorBehavior::forceSequential) or autopas_get_num_threads() == 1)) {
      cellIndex = 1;
    }
    do {
      // Increment the particle index. If this breaches the end of a cell, go to the next one and reset particleIndex.
      if (++particleIndex >= this->_cells[cellIndex].numParticles()) {
        cellIndex += stride;
        particleIndex = 0;
      }
      // If we notice that there is nothing else to look at set invalid values, so we get a nullptr next time and break.
      if (cellIndex > (iteratorBehavior & IteratorBehavior::owned) ? 0 : this->_cells.size() - 1) {
        cellIndex = std::numeric_limits<size_t>::max();
        break;
      }
      // Repeat this as long as the current particle is not interesting.
      //  - coordinates are in region of interest
      //  - ownership fits to the iterator behavior
    } while (not utils::inBox(this->_cells[cellIndex][particleIndex].getR(), boxMin, boxMax) or
             not(static_cast<unsigned int>(this->_cells[cellIndex][particleIndex].getOwnershipState()) &
                 static_cast<unsigned int>(iteratorBehavior)));

    return {retPtr, cellIndex, particleIndex};
  };

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

  ParticleCell &getCell() { return this->_cells[0]; };

  ParticleCell &getHaloCell() { return this->_cells[1]; };
};

}  // namespace autopas
