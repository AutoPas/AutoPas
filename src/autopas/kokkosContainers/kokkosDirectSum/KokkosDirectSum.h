/**
 * @file KokkosCellBasedParticleContainer.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include "autopas/kokkosContainers/KokkosCellBasedParticleContainer.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/kokkosContainers/kokkosDirectSum/traversals/KokkosDSTraversalInterface.h"
#include "autopas/kokkosContainers/KokkosCellPairTraversals/KokkosCellPairTraversal.h"

namespace autopas {

/**
 * KokkosDirectSum class.
 */
template <class Particle>
class KokkosDirectSum : public KokkosCellBasedParticleContainer<Particle> {
  using ParticleCell = KokkosParticleCell<Particle>;

 public:
  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  KokkosDirectSum(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin)
      : KokkosCellBasedParticleContainer<Particle>(boxMin, boxMax, cutoff, skin, 2) {}

  ~KokkosDirectSum() = default;

  ContainerOption getContainerType() const override { return ContainerOption::kokkosDirectSum; };

  CellType getParticleCellTypeEnum() override { return CellType::KokkosCell; }

  bool getIsDirty() { return this->_isDirty; }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const Particle &p) override {
    this->_particles.addParticle(p);
    // TODO lgaertner: maybe check for existence of halo particles?
    this->_isDirty = true;
  }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addHaloParticleImpl(const Particle &p) override {
    Particle pCopy = p;
    pCopy.setOwnershipState(OwnershipState::halo);
    this->_particles.addHaloParticle(pCopy);
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
    auto *traversalInterface = dynamic_cast<KokkosDSTraversalInterface<ParticleCell> *>(traversal);
    auto *cellPairTraversal = dynamic_cast<KokkosCellPairTraversal<ParticleCell> *>(traversal);

    if (traversalInterface && cellPairTraversal) {
      cellPairTraversal->setCellsToTraverse(this->_cells);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in KokkosDirectSum::iteratePairwise");
    }

    cellPairTraversal->initTraversal();
    cellPairTraversal->traverseParticlePairs();
    cellPairTraversal->endTraversal();

  }

  std::vector<Particle> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      //      basically do nothing and return
    }

    //    TODO lgaertner
    std::vector<Particle> invalidParticles{};
//    this->_particles.forEach();
//    for (auto iter = getCell().begin(); iter.isValid(); ++iter) {
//      if (utils::notInBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
//        invalidParticles.push_back(*iter);
//        internal::deleteParticle(iter);
//      }
//    }

    this->deleteHaloParticles();
    this->_particles.binParticles([&](Particle &p) -> size_t { return p.isOwned() ? _OWNED : _HALO; }, this->_cells,
                                  "KokkosDirectSum::updateContainer:");
    this->_isDirty = false;

    return invalidParticles;
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
        this->_particles.forEach(forEachLambda, this->_cells[_OWNED], "KokkosDirectSum:forEach");
      } else if (behavior & IteratorBehavior::halo) {
        this->_particles.forEach(forEachLambda, this->_cells[_HALO], "KokkosDirectSum:forEach");
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
          this->_particles.reduce(reduceLambda, result, this->_cells[_OWNED], "KokkosDirectSum:reduce");
        } else if (behavior & IteratorBehavior::halo) {
          this->_particles.reduce(reduceLambda, result, this->_cells[_HALO], "KokkosDirectSum:reduce");
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
  size_t _OWNED{0};
  size_t _HALO{1};
};

}  // namespace autopas