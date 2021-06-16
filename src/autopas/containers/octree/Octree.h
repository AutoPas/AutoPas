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
   */
  Octree(std::array<double, 3> boxMin, std::array<double, 3> boxMax, const double cutoff, const double skin)
      : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skin) {
    // TODO(johannes): Obtain this from a configuration
    int unsigned treeSplitThreshold = 16;
    this->_cells.push_back(OctreeNodeWrapper<Particle>(boxMin, boxMax, treeSplitThreshold));

    // Extend the halo region with cutoff + skin in all dimensions
    auto interactionLength = cutoff + skin;
    auto haloBoxMin = utils::ArrayMath::subScalar(boxMin, interactionLength);
    auto haloBoxMax = utils::ArrayMath::addScalar(boxMax, interactionLength);
    this->_cells.push_back(OctreeNodeWrapper<Particle>(haloBoxMin, haloBoxMax, treeSplitThreshold));
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer() override {
    // This is a very primitive and inefficient way to recreate the container:
    // 1. Copy all particles out of the container
    // 2. Clear the container
    // 3. Insert the particles back into the container

    // leaving: all outside boxMin/Max

    std::vector<Particle> particles;
    this->_cells[CellTypes::OWNED].appendAllParticles(particles);

    deleteAllParticles();

    for (auto &particle : particles) {
      addParticleImpl(particle);
    }

    // logger.logTree(_root);

    auto result = std::vector<ParticleType>();
    return result;
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    printf("Johannes' Octree::iteratePairwise\n");

    auto *traversalInterface = dynamic_cast<OTTraversalInterface<ParticleCell> *>(traversal);
    if (traversalInterface) {
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
    this->_cells[CellTypes::HALO].addParticle(haloParticle);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    printf("Johannes' Octree::updateHaloParticle\n");
    return true;
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    printf("Johannes' Octree::rebuildNeighborLists\n");
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> begin(IteratorBehavior behavior) override {
    printf("Johannes' Octree::begin<..., true>\n");

    /*std::vector<OctreeLeafNode<Particle> *> leaves;
    _root->appendAllLeafNodesInside(leaves, this->getBoxMin(), this->getBoxMax());

    std::vector<OctreeLeafNode<Particle>> flatLeaves;
    for(auto *leaf : leaves) {
        flatLeaves.push_back(*leaf);
    }*/
    // std::vector<FullParticleCell<Particle>> cells;
    //_root->appendAllParticleCellsInside(cells);
    /*for(auto *leaf : leaves) {
        auto *cell = dynamic_cast<FullParticleCell<Particle> *>(leaf);
        cells.push_back(*cell);
    }*/
    // std::vector<FullParticleCell<Particle>> cells = flatLeaves;

    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::ParticleIterator<ParticleType, ParticleCell, true>(&this->_cells, 0, this, behavior));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> begin(IteratorBehavior behavior) const override {
    // printf("Johannes' Octree::begin<..., false>\n");
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::ParticleIterator<ParticleType, ParticleCell, false>(&this->_cells, 0, this, behavior));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                              const std::array<double, 3> &higherCorner,
                                                                              IteratorBehavior behavior) override {
    printf("Johannes' Octree::getRegionIterator<..., true>\n");
    // TODO(johannes): This is a bad implementation, it does not utilize the spacial structure of the octree :(
    std::vector<size_t> cellsOfInterest = {CellTypes::OWNED};  // TODO(johannes): Add the second cell here
    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, true>(&this->_cells, lowerCorner, higherCorner,
                                                                               cellsOfInterest, this, behavior));
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) const override {
    printf("Johannes' Octree::getRegionIterator<..., false>\n");
    // TODO(johannes): This is a bad implementation, it does not utilize the spacial structure of the octree :(
    std::vector<size_t> cellsOfInterest = {CellTypes::OWNED};  // TODO(johannes): Add the second cell here
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, false>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, this, behavior));
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // TODO(johannes): Figure out what these values should be
    std::array<unsigned long, 3> dims = {1, 1, 1};
    std::array<double, 3> cellLength = {1, 1, 1};
    return TraversalSelectorInfo(dims, 0.0, cellLength, 1);
  }

  /**
   * Get the number of particles that belong to this octree. (Only the owned particles.)
   *
   * @return The integer # of particles in the container
   */
  [[nodiscard]] unsigned long getNumParticles() const override { return this->_cells[CellTypes::OWNED].numParticles(); }

  /**
   * Deletes all particles from the container.
   */
  void deleteAllParticles() override { this->_cells[CellTypes::OWNED].clear(); }

  void deleteHaloParticles() override { this->_cells[CellTypes::HALO].clear(); }

  [[nodiscard]] bool cellCanContainHaloParticles(std::size_t i) const override {
    if (i > 1) {
      throw std::runtime_error("This cell container (octree) contains only two cells");
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
  OctreeLogger logger;
};

}  // namespace autopas
