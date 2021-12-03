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
  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  KokkosDirectSum(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin)
      : KokkosCellBasedParticleContainer<Particle>(boxMin, boxMax, cutoff, skin),
        _begin("KokkosDirectSum::_particles.begin", 2),
        _cellSize("KokkosDirectSum::_particles.cellSize", 2) {}

  ~KokkosDirectSum() = default;

  ContainerOption getContainerType() const override { return ContainerOption::kokkosDirectSum; };

  CellType getParticleCellTypeEnum() override { return CellType::KokkosCell; }

  bool getIsDirty() {
    return this->_isDirty;
  }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const Particle &p) override {
    this->_particles.addParticle(p);
    //TODO lgaertner: maybe check for existence of halo particles?
      this->_isDirty = true;
  }

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
    this->_particles.forEach([&](Particle &p) { p.setOwnershipState(OwnershipState::dummy); }, IteratorBehavior::halo);
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    // TODO lgaertner
  }

  std::vector<Particle> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      //      basically do nothing and return
    }

    this->deleteHaloParticles();
    this->_particles.binParticles([&](Particle &p) -> size_t { return p.isOwned() ? _OWNED : _HALO; }, _begin,
                                  _cellSize, "KokkosDirectSum::updateContainer: ");
    this->_isDirty = false;

    return std::vector<Particle>();
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    if (behavior == IteratorBehavior::ownedOrHalo) {
      this->_particles.forEach(forEachLambda, "KokkosDirectSum:forEach");
    } else {
      if (not this->_isDirty) {
        this->_particles.forEach(forEachLambda, behavior, "KokkosDirectSum:forEach(behavior)");
      }
      if (behavior & IteratorBehavior::owned) {
        this->_particles.forEach(forEachLambda, _begin[_OWNED], _cellSize[_OWNED], "KokkosDirectSum:forEach");
      } else if (behavior & IteratorBehavior::halo) {
        this->_particles.forEach(forEachLambda, _begin[_HALO], _cellSize[_HALO], "KokkosDirectSum:forEach");
      }
    }
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior) {
    if (behavior == IteratorBehavior::ownedOrHalo) {
      this->_particles.reduce(reduceLambda, result, "KokkosDirectSum:reduce");
    } else {
      if (this->_isDirty) {
        this->_particles.reduce(reduceLambda, result, behavior, "KokkosDirectSum:reduce(behavior)");
      } else {
        if (behavior & IteratorBehavior::owned) {
          this->_particles.reduce(reduceLambda, result, _begin[_OWNED], _cellSize[_OWNED], "KokkosDirectSum:reduce");
        } else if (behavior & IteratorBehavior::halo) {
          this->_particles.reduce(reduceLambda, result, _begin[_HALO], _cellSize[_HALO], "KokkosDirectSum:reduce");
        }
      }
    }
  }

  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    this->_particles.forEach(forEachLambda, lowerCorner, higherCorner, behavior,
                             "KokkosDirectSum:forEach(lowerCorner, higherCorner, behavior)");
  }

  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    this->_particles.reduce(reduceLambda, result, behavior, lowerCorner, higherCorner,
                            "KokkosDirectSum::reduce(behavior, lowerCorner, higherCorner)");
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
  Kokkos::View<size_t *> _begin;
  Kokkos::View<size_t *> _cellSize;

  size_t _OWNED{0};
  size_t _HALO{1};
};

}  // namespace autopas