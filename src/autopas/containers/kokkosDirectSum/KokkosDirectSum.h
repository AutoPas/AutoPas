/**
 * @file KokkosCellBasedParticleContainer.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include "autopas/containers/KokkosCellBasedParticleContainer.h"

namespace autopas {

/**
 * KokkosDirectSum class.
 */
template <class Particle>
class KokkosDirectSum : public KokkosCellBasedParticleContainer<Particle> {
 public:
  CellType getParticleCellTypeEnum() override { return CellType::KokkosCell; }

  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  KokkosDirectSum(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin)
      : KokkosCellBasedParticleContainer<Particle>(boxMin, boxMax, cutoff, skin) {}

  ~KokkosDirectSum() = default;

  ContainerOption getContainerType() const override { return ContainerOption::kokkosDirectSum; };

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const Particle &p) override { this->_particles.addParticle(p); }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addHaloParticleImpl(const Particle &p) override {
    Particle pCopy = p;
    pCopy.setOwnershipState(OwnershipState::halo);
    this->_particles.addParticle(pCopy);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle &haloParticle) override {
    // TODO lgaertner
    return false;
  }

  /**
   * @copydoc ParticleContainerInterface::rebuildNeighborLists()
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // TODO lgaertner
  }

  /**
   * @copydoc ParticleContainerInterface::deleteHaloParticles()
   */
  void deleteHaloParticles() override {
    // TODO lgaertner
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    // TODO lgaertner
  }

  std::vector<Particle> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      //      basically do nothing and return
    }

    return std::vector<Particle>();
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    if (behavior == IteratorBehavior::ownedOrHalo) {
      this->_particles.forEach(forEachLambda);
    } else {
      //      if (not this->_isDirty) {
      //        this->updateContainer(false);
      //      }
      //      if (behavior & IteratorBehavior::owned) {
      //        this->_particles.forEach(forEachLambda, _start[_OWNED], _cellSize[_OWNED]);
      //      } else if (behavior & IteratorBehavior::halo) {
      //        this->_particles.forEach(forEachLambda, _start[_HALO], _cellSize[_HALO]);
      //      }
      this->_particles.forEach(forEachLambda, behavior, "KokkosDirectSum:forEach");
    }
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior) {
    //    if (behavior == IteratorBehavior::ownedOrHalo) {
    //      this->_particles.reduce(reduceLambda, result);
    //    } else {
    //      if (behavior & IteratorBehavior::owned) {
    //        this->_particles.reduce(reduceLambda, result, _start[_OWNED], _cellSize[_OWNED]);
    //      } else if (behavior & IteratorBehavior::halo) {
    //        this->_particles.reduce(reduceLambda, result, _start[_HALO], _cellSize[_HALO]);
    //      }
    //    }
    this->_particles.reduce(reduceLambda, result, behavior, "KokkosDirectSum:reduce");
  }

  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
  }

  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    this->_particles.reduce(reduceLambda, result, behavior, lowerCorner, higherCorner, "KokkosDirectSum::reduce(behavior, lowerCorner, higherCorner)");
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // TODO lgaertner
    std::array<size_t, 3> dummy{};
    std::array<double, 3> dummy2{};

    return TraversalSelectorInfo(dummy, this->getInteractionLength(), dummy2, 0);
  }

 private:
  Kokkos::View<size_t[2]> _start;
  Kokkos::View<size_t[2]> _cellSize;

  size_t _OWNED{0};
  size_t _HALO{1};
};

}  // namespace autopas