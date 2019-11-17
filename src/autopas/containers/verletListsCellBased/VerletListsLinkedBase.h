/**
 * @file VerletListsLinkedBase.h
 * @author nguyen
 * @date 17.12.18
 */

#pragma once

#include "autopas/containers/ParticleContainer.h"
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
 * @tparam LinkedSoAArraysType SoAArraysType used by the linked cells container
 */
template <class Particle, class LinkedParticleCell, class LinkedSoAArraysType = typename Particle::SoAArraysType>
class VerletListsLinkedBase : public ParticleContainerInterface<FullParticleCell<Particle>> {
  using ParticleCell = FullParticleCell<Particle>;

 public:
  /**
   * Constructor of the VerletListsLinkedBase class.
   * The neighbor lists are build using a search radius of cutoff + skin.
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
   * @copydoc autopas::ParticleContainerInterface::addParticle
   * @note This function invalidates the neighbor lists.
   */
  void addParticle(Particle &p) override {
    _neighborListIsValid = false;
    _linkedCells.addParticle(p);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::addHaloParticle
   * @note This function invalidates the neighbor lists.
   */
  void addHaloParticle(Particle &haloParticle) override {
    _neighborListIsValid = false;
    _linkedCells.addHaloParticle(haloParticle);
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
  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<Particle> updateContainer() override {
    AutoPasLog(debug, "updating container");
    _neighborListIsValid = false;
    return _linkedCells.updateContainer();
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::isContainerUpdateNeeded()
   */
  bool isContainerUpdateNeeded() const override {
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // TODO: find a sensible value for ???
#pragma omp parallel for shared(outlierFound)  // if (this->_cells.size() / omp_get_max_threads() > ???)
#endif
    for (size_t cellIndex1d = 0; cellIndex1d < _linkedCells.getCells().size(); ++cellIndex1d) {
      std::array<double, 3> boxmin{0., 0., 0.};
      std::array<double, 3> boxmax{0., 0., 0.};
      _linkedCells.getCellBlock().getCellBoundingBox(cellIndex1d, boxmin, boxmax);
      boxmin = utils::ArrayMath::subScalar(boxmin, this->getSkin() / 2.);
      boxmax = utils::ArrayMath::addScalar(boxmax, this->getSkin() / 2.);
      for (auto iter = _linkedCells.getCells()[cellIndex1d].begin(); iter.isValid(); ++iter) {
        if (not utils::inBox(iter->getR(), boxmin, boxmax)) {
          outlierFound = true;  // we need an update
          break;
        }
      }
      if (outlierFound) cellIndex1d = _linkedCells.getCells().size();
    }
    if (outlierFound) {
      AutoPasLog(debug,
                 "VerletLists: containerUpdate needed! Particles are fast. You "
                 "might want to increase the skin radius or decrease the rebuild "
                 "frequency.");
    } else {
      AutoPasLog(debug,
                 "VerletLists: containerUpdate not yet needed. Particles are slow "
                 "enough.");
    }
    return outlierFound;
  }

  /**
   * Searches the provided halo particle and updates the found particle.
   * Searches for the provided particle within the halo cells of the container
   * and overwrites the found particle with the provided particle.
   * @param particle
   * @return true if a particle was found and updated, false if it was not found.
   */
  bool updateHaloParticle(Particle &particle) override {
    Particle pCopy = particle;
    pCopy.setOwned(false);
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
  ParticleIteratorWrapper<Particle, true> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return _linkedCells.begin(behavior);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  ParticleIteratorWrapper<Particle, false> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    return _linkedCells.begin(behavior);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  ParticleIteratorWrapper<Particle, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return _linkedCells.getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  ParticleIteratorWrapper<Particle, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    return _linkedCells.getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * Get the dimension of the used cellblock including the haloboxes.
   * @return the dimensions of the used cellblock
   */
  const std::array<std::size_t, 3> &getCellsPerDimension() const {
    return _linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
  }

  /**
   * Generates a traversal selector info for this container.
   * @return Traversal selector info for this container.
   */
  TraversalSelectorInfo getTraversalSelectorInfo() const override {
    return TraversalSelectorInfo(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                 this->getInteractionLength(), this->_linkedCells.getCellBlock().getCellLength());
  }

  [[nodiscard]] std::unique_ptr<fmm::FmmTree> getFastMultipoleMethodTree() override {
    return std::move(_linkedCells.getFastMultipoleMethodTree());
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMax()
   */
  const std::array<double, 3> &getBoxMax() const override final { return _linkedCells.getBoxMax(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMax()
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override final { _linkedCells.setBoxMax(boxMax); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMin()
   */
  const std::array<double, 3> &getBoxMin() const override final { return _linkedCells.getBoxMin(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMin()
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override final { _linkedCells.setBoxMin(boxMin); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getCutoff()
   */
  double getCutoff() const override final { return _linkedCells.getCutoff(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) override final { _linkedCells.setCutoff(cutoff); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getSkin()
   */
  double getSkin() const override final { return _linkedCells.getSkin(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setSkin()
   */
  void setSkin(double skin) override final { _linkedCells.setSkin(skin); }

  /**
   * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
   */
  double getInteractionLength() const override final { return _linkedCells.getInteractionLength(); }

 protected:
  /// internal linked cells storage, handles Particle storage and used to build verlet lists
  LinkedCells<LinkedParticleCell, LinkedSoAArraysType> _linkedCells;

  /// specifies if the neighbor list is currently valid
  bool _neighborListIsValid{false};

  /// specifies if the current verlet list was built for newton3
  bool _verletBuiltNewton3{false};
};

}  // namespace autopas
