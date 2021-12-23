/**
 * @file VerletListsLinkedBase.h
 * @author nguyen
 * @date 17.12.18
 */

#pragma once

#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ParticleCellHelpers.h"

namespace autopas {

/**
 * Base class for Verlet lists which use an underlying linked cells container.
 * Implementation have to use a constant cutoff radius of the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle
 * @tparam LinkedParticleCells ParticleCells used by the linked cells container
 */
template <class Particle>
class VerletListsLinkedBase : public ParticleContainerInterface<Particle> {
 public:
  /**
   * Constructor of the VerletListsLinkedBase class.
   * The neighbor lists are build using a search radius of cutoff + skin.LinkedParticleCell::ParticleType
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param applicableTraversals all applicable traversals
   * @param cellSizeFactor cell size factor relative to cutoff. Verlet lists are only implemented for values >= 1.0
   * (smaller values are set to 1.0).
   */
  VerletListsLinkedBase(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                        const double skin, const std::set<TraversalOption> &applicableTraversals,
                        const double cellSizeFactor)
      : _linkedCells(boxMin, boxMax, cutoff, skin, std::max(1.0, cellSizeFactor)) {
    if (cellSizeFactor < 1.0) {
      AutoPasLog(debug, "VerletListsLinkedBase: CellSizeFactor smaller 1 detected. Set to 1.");
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getParticleCellTypeEnum()
   */
  CellType getParticleCellTypeEnum() override { return _linkedCells.getParticleCellTypeEnum(); };

  /**
   * @copydoc autopas::ParticleContainerInterface::addParticleImpl
   * @note This function invalidates the neighbor lists.
   */
  void addParticleImpl(const Particle &p) override {
    _neighborListIsValid = false;
    // position is already checked, so call impl directly.
    _linkedCells.addParticleImpl(p);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::addHaloParticleImpl
   * @note This function invalidates the neighbor lists.
   */
  void addHaloParticleImpl(const Particle &haloParticle) override {
    _neighborListIsValid = false;
    // position is already checked, so call impl directly.
    _linkedCells.addHaloParticleImpl(haloParticle);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getNumParticles()
   */
  unsigned long getNumParticles() const override { return _linkedCells.getNumParticles(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::deleteHaloParticles
   * @note This function invalidates the neighbor lists.
   */
  void deleteHaloParticles() override {
    _neighborListIsValid = false;
    _linkedCells.deleteHaloParticles();
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::deleteAllParticles
   * @note This function invalidates the neighbor lists.
   */
  void deleteAllParticles() override {
    _neighborListIsValid = false;
    _linkedCells.deleteAllParticles();
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::updateContainer()
   * @note This function invalidates the neighbor lists.
   */
  [[nodiscard]] std::vector<Particle> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return autopas::LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(_linkedCells);
    }
    _neighborListIsValid = false;
    return _linkedCells.updateContainer(false);
  }

  /**
   * Searches the provided halo particle and updates the found particle.
   * Searches for the provided particle within the halo cells of the container
   * and overwrites the found particle with the provided particle.
   * @param particle
   * @return true if a particle was found and updated, false if it was not found.
   */
  bool updateHaloParticle(const Particle &particle) override {
    Particle pCopy = particle;
    pCopy.setOwnershipState(OwnershipState::halo);
    auto cells = _linkedCells.getCellBlock().getNearbyHaloCells(pCopy.getR(), this->getSkin());
    for (auto cellptr : cells) {
      bool updated = internal::checkParticleInCellAndUpdateByID(*cellptr, pCopy);
      if (updated) {
        return true;
      }
    }
    AutoPasLog(trace,
               "updateHaloParticle was not able to update particle at "
               "[{}, {}, {}]",
               pCopy.getR()[0], pCopy.getR()[1], pCopy.getR()[2]);
    return false;
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, true> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) override {
    return _linkedCells.begin(behavior);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const override {
    return _linkedCells.begin(behavior);
  }

  /**
   * @copydoc autopas::LinkedCells::forEach()
   */
  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    return _linkedCells.template forEach<parallel>(forEachLambda, behavior);
  }

  /**
   * @copydoc autopas::LinkedCells::reduce()
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior) {
    return _linkedCells.reduce(reduceLambda, result, behavior);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                          const std::array<double, 3> &higherCorner,
                                                                          IteratorBehavior behavior) override {
    return _linkedCells.getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, false> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                           const std::array<double, 3> &higherCorner,
                                                                           IteratorBehavior behavior) const override {
    return _linkedCells.getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc autopas::LinkedCells::forEachInRegion()
   */
  template <bool parallel, typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    _linkedCells.template forEachInRegion<parallel>(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc autopas::LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    _linkedCells.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
  }

  /**
   * Get the dimension of the used cellblock including the haloboxes.
   * @return the dimensions of the used cellblock
   */
  [[nodiscard]] const std::array<std::size_t, 3> &getCellsPerDimension() const {
    return _linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
  }

  /**
   * Generates a traversal selector info for this container.
   * @return Traversal selector info for this container.
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    return TraversalSelectorInfo(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                 this->getInteractionLength(), this->_linkedCells.getCellBlock().getCellLength(), 0);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMax()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override final { return _linkedCells.getBoxMax(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMax()
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override final { _linkedCells.setBoxMax(boxMax); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMin()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override final { return _linkedCells.getBoxMin(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMin()
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override final { _linkedCells.setBoxMin(boxMin); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getCutoff()
   */
  [[nodiscard]] double getCutoff() const override final { return _linkedCells.getCutoff(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) override final { _linkedCells.setCutoff(cutoff); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getSkin()
   */
  [[nodiscard]] double getSkin() const override final { return _linkedCells.getSkin(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setSkin()
   */
  void setSkin(double skin) override final { _linkedCells.setSkin(skin); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
   */
  [[nodiscard]] double getInteractionLength() const override final { return _linkedCells.getInteractionLength(); }

 protected:
  /// internal linked cells storage, handles Particle storage and used to build verlet lists
  LinkedCells<Particle> _linkedCells;

  /// specifies if the neighbor list is currently valid
  bool _neighborListIsValid{false};

  /// specifies if the current verlet list was built for newton3
  bool _verletBuiltNewton3{false};
};

}  // namespace autopas
