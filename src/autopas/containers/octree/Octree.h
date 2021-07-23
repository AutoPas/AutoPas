/**
 * @file Octree.h
 *
 * @author Johannes Spies
 * @date 09.04.2021
 */

#pragma once

#include <cstdio>
#include <list>
#include <stack>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/octree/OctreeLeafNode.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeNodeWrapper.h"
#include "autopas/containers/octree/traversals/OTTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/utils/ParticleCellHelpers.h"
#include "autopas/utils/logging/OctreeLogger.h"

namespace autopas {

/**
 * The octree is a CellBasedParticleContainer that is comprised of two cells:
 *
 * @tparam Particle
 */
template <class Particle>
class Octree : public CellBasedParticleContainer<OctreeNodeWrapper<Particle>>,
               public internal::CellBorderAndFlagManager {
 public:
  /**
   * The particle cell used in this CellBasedParticleContainer
   */
  using ParticleCell = OctreeNodeWrapper<Particle>;

  /**
   * The particle type used in this container.
   */
  using ParticleType = typename ParticleCell::ParticleType;

  /**
   * This particle container contains two cells. Both cells are octrees.
   * - this->_cells[CellTypes::OWNED] yields the octree that contains all owned particles
   * - this->_cells[CellTypes::HALO] yields the octree that contains all halo particles
   */
  enum CellTypes { OWNED = 0, HALO = 1 };

  /**
   * Construct a new octree with two sub-octrees: One for the owned particles and one for the halo particles.
   * @param boxMin The minimum coordinate of the enclosing box
   * @param boxMax The maximum coordinate of the enclosing box
   * @param cutoff The cutoff radius
   * @param skin The skin radius
   * @param cellSizeFactor The cell size factor
   */
  Octree(std::array<double, 3> boxMin, std::array<double, 3> boxMax, const double cutoff, const double skin,
         const double cellSizeFactor)
      : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skin) {
    // @todo Obtain this from a configuration, reported in https://github.com/AutoPas/AutoPas/issues/624
    int unsigned treeSplitThreshold = 16;

    double interactionLength = this->getInteractionLength();

    // Create the octree for the owned particles
    this->_cells.push_back(
        OctreeNodeWrapper<Particle>(boxMin, boxMax, treeSplitThreshold, interactionLength, cellSizeFactor));

    // Extend the halo region with cutoff + skin in all dimensions
    auto haloBoxMin = utils::ArrayMath::subScalar(boxMin, interactionLength);
    auto haloBoxMax = utils::ArrayMath::addScalar(boxMax, interactionLength);
    // Create the octree for the halo particles
    this->_cells.push_back(
        OctreeNodeWrapper<Particle>(haloBoxMin, haloBoxMax, treeSplitThreshold, interactionLength, cellSizeFactor));
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer() override {
    // This is a very primitive and inefficient way to rebuild the container:
    // 1. Copy all particles out of the container
    // 2. Clear the container
    // 3. Insert the particles back into the container

    // leaving: all outside boxMin/Max

    // @todo Make this less indirect. (Find a better way to iterate all particles inside the octree to change
    //   this function back to a function that actually copies all particles out of the octree.)
    //   The problem is captured by https://github.com/AutoPas/AutoPas/issues/622
    std::vector<Particle *> particleRefs;
    this->_cells[CellTypes::OWNED].appendAllParticles(particleRefs);
    std::vector<Particle> particles{};
    particles.reserve(particleRefs.size());
    std::vector<Particle> invalidParticles{};
    for (auto *p : particleRefs) {
      if (utils::inBox(p->getR(), this->getBoxMin(), this->getBoxMax())) {
        particles.push_back(*p);
      } else {
        invalidParticles.push_back(*p);
      }
    }
    this->deleteAllParticles();

    for (auto &particle : particles) {
      addParticleImpl(particle);
    }

    return invalidParticles;
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    if (auto *traversalInterface = dynamic_cast<OTTraversalInterface<ParticleCell> *>(traversal)) {
      traversalInterface->setCells(&this->_cells);
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::octree; }

  /**
   * @copydoc ParticleContainerInterface::getParticleCellTypeEnum()
   */
  [[nodiscard]] CellType getParticleCellTypeEnum() override { return CellType::FullParticleCell; }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const ParticleType &p) override { this->_cells[CellTypes::OWNED].addParticle(p); }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    ParticleType p_copy = haloParticle;
    p_copy.setOwnershipState(OwnershipState::halo);
    this->_cells[CellTypes::HALO].addParticle(p_copy);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwnershipState(OwnershipState::halo);
    return internal::checkParticleInCellAndUpdateByIDAndPosition(this->_cells[CellTypes::HALO], pCopy, this->getSkin());
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {}

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> begin(IteratorBehavior behavior) override {
    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::ParticleIterator<ParticleType, ParticleCell, true>(&this->_cells, 0, this, behavior));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> begin(IteratorBehavior behavior) const override {
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::ParticleIterator<ParticleType, ParticleCell, false>(&this->_cells, 0, this, behavior));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                              const std::array<double, 3> &higherCorner,
                                                                              IteratorBehavior behavior) override {
    std::vector<size_t> cellsOfInterest;
    if (behavior & IteratorBehavior::owned) {
      cellsOfInterest.push_back(CellTypes::OWNED);
    }
    if (behavior & IteratorBehavior::halo) {
      cellsOfInterest.push_back(CellTypes::HALO);
    }
    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, true>(&this->_cells, lowerCorner, higherCorner,
                                                                               cellsOfInterest, this, behavior));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) const override {
    std::vector<size_t> cellsOfInterest;
    if (behavior & IteratorBehavior::owned) {
      cellsOfInterest.push_back(CellTypes::OWNED);
    }
    if (behavior & IteratorBehavior::halo) {
      cellsOfInterest.push_back(CellTypes::HALO);
    }
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, false>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, this, behavior));
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // this is a dummy since it is not actually used
    std::array<unsigned long, 3> dims = {1, 1, 1};
    std::array<double, 3> cellLength = utils::ArrayMath::sub(this->getBoxMax(), this->getBoxMin());
    return TraversalSelectorInfo(dims, this->getInteractionLength(), cellLength, 0);
  }

  /**
   * Get the number of particles that belong to this octree. (Owned and and halo.)
   *
   * @return The integer # of particles in the container
   */
  [[nodiscard]] unsigned long getNumParticles() const override {
    return this->_cells[CellTypes::OWNED].numParticles() + this->_cells[CellTypes::HALO].numParticles();
  }

  void deleteHaloParticles() override { this->_cells[CellTypes::HALO].clear(); }

  [[nodiscard]] bool cellCanContainHaloParticles(std::size_t i) const override {
    if (i > 1) {
      throw std::runtime_error("[Octree.h]: This cell container (octree) contains only two cells");
    }
    return i == CellTypes::HALO;
  }

  [[nodiscard]] bool cellCanContainOwnedParticles(std::size_t i) const override {
    if (i > 1) {
      throw std::runtime_error("[Octree.h]: This cell container (octree) contains only two cells");
    }
    return i == CellTypes::OWNED;
  }

 private:
  /**
   * A logger that can be called to log the octree data structure.
   */
  OctreeLogger<Particle> logger;
};

}  // namespace autopas
