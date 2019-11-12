/**
 * @file DirectSum.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/directSum/DirectSumTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/AutoPasMacros.h"
#include "autopas/utils/CudaStreamHandler.h"
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
template <class ParticleCell>
class DirectSum : public ParticleContainer<ParticleCell> {
 public:
  /**
   *  Type of the Particle.
   */
  using ParticleType = typename ParticleContainer<ParticleCell>::ParticleType;

  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  DirectSum(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin)
      : ParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skin), _cellBorderFlagManager() {
    this->_cells.resize(2);
  }

  ContainerOption getContainerType() const override { return ContainerOption::directSum; }

  /**
   * @copydoc ParticleContainerInterface::addParticle()
   */
  void addParticle(ParticleType &p) override {
    if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      getCell().addParticle(p);
    } else {
      utils::ExceptionHandler::exception("DirectSum: trying to add a particle that is not in the bounding box.\n" +
                                         p.toString());
    }
  }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticle()
   */
  void addHaloParticle(ParticleType &haloParticle) override {
    ParticleType p_copy = haloParticle;
    p_copy.setOwned(false);
    getHaloCell().addParticle(p_copy);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwned(false);
    return internal::checkParticleInCellAndUpdateByIDAndPosition(getHaloCell(), pCopy, this->getSkin());
  }

  void deleteHaloParticles() override { getHaloCell().clear(); }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    AutoPasLog(debug, "Using traversal {}.", traversal->getTraversalType().to_string());

    // Check if traversal is allowed for this container and give it the data it needs.
    auto *traversalInterface = dynamic_cast<DirectSumTraversalInterface<ParticleCell> *>(traversal);
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

  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<ParticleType> updateContainer() override {
    // first we delete halo particles, as we don't want them here.
    deleteHaloParticles();
    std::vector<ParticleType> invalidParticles{};
    for (auto iter = getCell().begin(); iter.isValid(); ++iter) {
      if (utils::notInBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        invalidParticles.push_back(*iter);
        iter.deleteCurrentParticle();
      }
    }
    return invalidParticles;
  }

  bool isContainerUpdateNeeded() const override {
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // @todo: find a sensible value for ???
#pragma omp parallel shared(outlierFound)  // if (this->_cells.size() / omp_get_max_threads() > ???)
#endif
    for (auto iter = this->begin(); iter.isValid() && (not outlierFound); ++iter) {
      if (utils::notInBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        outlierFound = true;
      }
    }
    return outlierFound;
  }

  TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // direct sum technically consists of two cells (owned + halo)
    return TraversalSelectorInfo(
        {2, 0, 0},
        this->getCutoff() /*intentionally use cutoff here, as the directsumtraversal should be using the cutoff.*/,
        utils::ArrayMath::sub(this->getBoxMax(), this->getBoxMin()));
  }

  ParticleIteratorWrapper<ParticleType, true> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<ParticleType, true>(new internal::ParticleIterator<ParticleType, ParticleCell, true>(
        &this->_cells, 0, &_cellBorderFlagManager, behavior));
  }

  ParticleIteratorWrapper<ParticleType, false> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::ParticleIterator<ParticleType, ParticleCell, false>(&this->_cells, 0, &_cellBorderFlagManager,
                                                                          behavior));
  }

  ParticleIteratorWrapper<ParticleType, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    std::vector<size_t> cellsOfInterest;

    switch (behavior) {
      case IteratorBehavior::ownedOnly:
        cellsOfInterest.push_back(0);
        break;
      case IteratorBehavior::haloOnly:
        // for haloOnly all cells can contain halo particles!
        [[fallthrough]];
      case IteratorBehavior::haloAndOwned:
        cellsOfInterest.push_back(0);
        cellsOfInterest.push_back(1);
        break;
    }

    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, true>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBorderFlagManager, behavior));
  }

  ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    std::vector<size_t> cellsOfInterest;

    switch (behavior) {
      case IteratorBehavior::ownedOnly:
        cellsOfInterest.push_back(0);
        break;
      case IteratorBehavior::haloOnly:
        // for haloOnly all cells can contain halo particles!
        [[fallthrough]];
      case IteratorBehavior::haloAndOwned:
        cellsOfInterest.push_back(0);
        cellsOfInterest.push_back(1);
        break;
    }

    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, false>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBorderFlagManager, behavior));
  }

 private:
  class DirectSumCellBorderAndFlagManager : public internal::CellBorderAndFlagManager {
    /**
     * the index type to access the particle cells
     */
    typedef std::size_t index_t;

   public:
    bool cellCanContainHaloParticles(index_t index1d) const override { return index1d == 1; }

    bool cellCanContainOwnedParticles(index_t index1d) const override { return index1d == 0; }

  } _cellBorderFlagManager;

  ParticleCell &getCell() { return this->_cells.at(0); };

  ParticleCell &getHaloCell() { return this->_cells.at(1); };
};

}  // namespace autopas
