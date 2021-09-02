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
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/directSum/traversals/DSTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
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
   * @param skin
   */
  DirectSum(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin)
      : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skin), _cellBorderFlagManager() {
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
    return internal::checkParticleInCellAndUpdateByIDAndPosition(getHaloCell(), pCopy, this->getSkin());
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

  [[nodiscard]] std::vector<ParticleType> updateContainer() override {
    // first we delete halo particles, as we don't want them here.
    deleteHaloParticles();
    getCell().deleteDummyParticles();

    std::vector<ParticleType> invalidParticles{};
    for (auto iter = getCell().begin(); iter.isValid(); ++iter) {
      if (utils::notInBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        invalidParticles.push_back(*iter);
        internal::deleteParticle(iter);
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

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) override {
    return ParticleIteratorWrapper<ParticleType, true>(new internal::ParticleIterator<ParticleType, ParticleCell, true>(
        &this->_cells, 0, &_cellBorderFlagManager, behavior, nullptr));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const override {
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::ParticleIterator<ParticleType, ParticleCell, false>(&this->_cells, 0, &_cellBorderFlagManager,
                                                                          behavior, nullptr));
  }

  /**
   * @copydoc LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    auto forEach = [&](ParticleCell cell) {
      for (Particle p : cell._particles) {
        forEachLambda(p);
      }
    };

    std::vector<size_t> cellsOfInterest;

    if (behavior & IteratorBehavior::owned) {
      getCell().forEach(forEachLambda);
      cellsOfInterest.push_back(0);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().forEach(forEachLambda);
      cellsOfInterest.push_back(1);
    }
    // sanity check
    if (cellsOfInterest.empty()) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  /**
   * @copydoc LinkedCells::reduce()
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior) {
    auto forEach = [&](ParticleCell cell) {
      for (Particle p : cell._particles) {
        reduceLambda(p, result);
      }
    };

    std::vector<size_t> cellsOfInterest;

    if (behavior & IteratorBehavior::owned) {
      getCell().reduce(reduceLambda, result);
      cellsOfInterest.push_back(0);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().reduce(reduceLambda, result);
      cellsOfInterest.push_back(1);
    }
    // sanity check
    if (cellsOfInterest.empty()) {
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
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBorderFlagManager, behavior, nullptr));
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
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBorderFlagManager, behavior, nullptr));
  }

  /**
   * @copydoc LinkedCells::forEachInRegion()
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    std::vector<size_t> cellsOfInterest;

    if (behavior & IteratorBehavior::owned) {
      getCell().forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      cellsOfInterest.push_back(0);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      cellsOfInterest.push_back(1);
    }

    // sanity check
    if (cellsOfInterest.empty()) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  /**
   * @copydoc LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reductionLambda, A &reductionValue, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    std::vector<size_t> cellsOfInterest;

    if (behavior & IteratorBehavior::owned) {
      getCell().reduce(reductionLambda, reductionValue, lowerCorner, higherCorner, behavior);
      cellsOfInterest.push_back(0);
    }
    if (behavior & IteratorBehavior::halo) {
      getHaloCell().reduce(reductionLambda, reductionValue, lowerCorner, higherCorner, behavior);
      cellsOfInterest.push_back(1);
    }

    // sanity check
    if (cellsOfInterest.empty()) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
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

  ParticleCell &getCell() { return this->_cells.at(0); };

  ParticleCell &getHaloCell() { return this->_cells.at(1); };
};

}  // namespace autopas
