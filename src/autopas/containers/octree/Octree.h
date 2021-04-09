/**
 * @file Octree.h
 *
 * @author Johannes Spies
 * @date 09.04.2021
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include <cstdio>

namespace autopas {

// TODO(johannes): Documentation
template <class Particle>
class Octree : public CellBasedParticleContainer<FullParticleCell<Particle>> {
 public:
  using ParticleCell = FullParticleCell<Particle>;
  using ParticleType = typename ParticleCell::ParticleType;

  Octree(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
         const double skin) : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skin) {
    printf("Johannes' Octree()\n");
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer() override {
    printf("Johannes' Octree::updateContainer\n");
    auto result = std::vector<ParticleType>();
    return result;
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    printf("Johannes' Octree::iteratePairwise\n");
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
  void addParticleImpl(const ParticleType &p) override {
    printf("Johannes' Octree::addParticleImpl\n");
  }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    printf("Johannes' Octree::addHaloParticleImpl\n");
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    printf("Johannes' Octree::updateHaloParticle\n");
    return true;
  }

  void deleteHaloParticles() override {
    printf("Johannes' Octree::deleteHaloParticles\n");
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    printf("Johannes' Octree::rebuildNeighborLists\n");
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> begin(IteratorBehavior behavior) override {
    printf("Johannes' Octree::begin<..., true>\n");
    return ParticleIteratorWrapper<ParticleType, true>();
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> begin(IteratorBehavior behavior) const override {
    printf("Johannes' Octree::begin<..., false>\n");
    return ParticleIteratorWrapper<ParticleType, false>();
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) override {
    printf("Johannes' Octree::getRegionIterator<..., true>\n");
    return ParticleIteratorWrapper<ParticleType, true>();
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) const override {
    printf("Johannes' Octree::getRegionIterator<..., false>\n");
    return ParticleIteratorWrapper<ParticleType, false>();
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
};

} // namespace autopas
