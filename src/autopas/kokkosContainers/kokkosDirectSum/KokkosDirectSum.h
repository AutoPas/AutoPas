/**
 * @file KokkosCellBasedParticleContainer.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include "autopas/cells/KokkosParticleCell.h"
#include "autopas/kokkosContainers/KokkosCellBasedParticleContainer.h"
#include "autopas/kokkosContainers/kokkosCellPairTraversals/KokkosCellPairTraversal.h"
#include "autopas/kokkosContainers/kokkosDirectSum/traversals/KokkosDSTraversalInterface.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
/**
 * KokkosDirectSum class.
 */
template <class Particle>
class KokkosDirectSum : public KokkosCellBasedParticleContainer<Particle> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = KokkosParticleCell<Particle>;

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

  CellType getParticleCellTypeEnum() override { return CellType::KokkosParticleCell; }

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
    this->_particles.addParticle(pCopy);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle &hp) override {
    Particle haloParticle = hp;
    haloParticle.setOwnershipState(OwnershipState::halo);

    bool isFound = false;
    auto lambda = [&](Particle &p) {
      if (p.isDummy() && p.getID() == haloParticle.getID()) {
        auto distanceVec = autopas::utils::ArrayMath::sub(p.getR(), haloParticle.getR());
        auto distanceSqr = autopas::utils::ArrayMath::dot(distanceVec, distanceVec);
        if (distanceSqr < this->getSkin() * this->getSkin()) {
          p = haloParticle;
          // found the particle, return true
          isFound = true;  // should not run into race conditioning problems
        }
      }
    };

    if (this->_isDirty) {
      this->_particles.template forEach<true>(lambda, IteratorBehavior::ownedOrHaloOrDummy,
                                              "KokkosDirectSum::updateHaloParticle");
    } else {
      this->_particles.template forEach<true>(lambda, this->_cells[_HALO]);
    }
    return isFound;
  }

  void inspectContainer() {
    auto view = this->getCellsHost();
    auto cellOwned = view[_OWNED];
    auto cellHalo = view[_HALO];
    this->_particles.inspect();
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
    // TODO update numParticles?
    this->_particles.template forEach<true>([&](Particle &p) { p.setOwnershipState(OwnershipState::dummy); },
                                            IteratorBehavior::halo);
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    //    inspectContainer();
    auto name = traversal->getTraversalType();
    auto cell = traversal->isApplicable();
    auto *traversalInterface = dynamic_cast<KokkosDSTraversalInterface<ParticleCell> *>(traversal);
    auto *cellPairTraversal = dynamic_cast<KokkosCellPairTraversal<ParticleCell> *>(traversal);
    if (traversalInterface && cellPairTraversal) {
      cellPairTraversal->setCellsToTraverse(this->_cells);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in KokkosDirectSum::iteratePairwise " +
          traversal->getTraversalType().to_string());
    }

    cellPairTraversal->initTraversal();
    cellPairTraversal->traverseParticlePairs();
    cellPairTraversal->endTraversal();
  }

  std::vector<Particle> updateContainer(bool keepNeighborListsValid) override {
    std::vector<Particle> invalidParticles{};

    //    TODO lgaertner
    this->deleteHaloParticles();

    auto boxMin = this->getBoxMin();
    auto boxMax = this->getBoxMax();

    this->_particles.template forEach<false>(
        [&](Particle &p) {
          if (utils::notInBox(p.getR(), boxMin, boxMax)) {
            invalidParticles.push_back(p);
            p.setOwnershipState(OwnershipState::dummy);
          }
        },
        IteratorBehavior::owned);

    if (not keepNeighborListsValid) {
      this->_particles.template binParticles<false>([&](Particle &p) -> size_t { return assignCellToParticle(p); },
                                                    this->_cells, "KokkosDirectSum::updateContainer:");
      this->_isDirty = false;
    }

    return invalidParticles;
  }

  void resortContainerAndDeleteDummies() {
    this->_particles.template binParticles<true>([&](Particle &p) -> size_t { return assignCellToParticle(p); },
                                                 this->_cells, "KokkosDirectSum::updateContainer:");

    inspectContainer();
    this->_isDirty = false;
  }

  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    if (behavior == IteratorBehavior::ownedOrHalo) {
      this->_particles.template forEach<parallel>(forEachLambda, "KokkosDirectSum:forEach");
    } else {
      if (this->_isDirty) {
        this->_particles.template forEach<parallel>(forEachLambda, behavior, "KokkosDirectSum:forEach(behavior)");
      }
      if (behavior & IteratorBehavior::owned) {
        this->_particles.template forEach<parallel>(forEachLambda, this->_cells[_OWNED], "KokkosDirectSum:forEach");
      } else if (behavior & IteratorBehavior::halo) {
        this->_particles.template forEach<parallel>(forEachLambda, this->_cells[_HALO], "KokkosDirectSum:forEach");
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

  template <bool parallel, typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    this->_particles.template forEach<parallel>(forEachLambda, lowerCorner, higherCorner, behavior,
                                                "KokkosDirectSum:forEach(lowerCorner, higherCorner, behavior)");
  }

  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    this->_particles.reduce(reduceLambda, result, behavior, lowerCorner, higherCorner,
                            "KokkosDirectSum::reduce(behavior, lowerCorner, higherCorner)");
  }

  size_t assignCellToParticle(Particle &p) override {
    return p.isOwned() ? _OWNED : _HALO;
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