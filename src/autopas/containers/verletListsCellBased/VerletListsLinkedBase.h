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
class VerletListsLinkedBase : public ParticleContainer<Particle, FullParticleCell<Particle>> {
  typedef FullParticleCell<Particle> ParticleCell;

 public:
  /**
   * Constructor of the VerletListsLinkedBase class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * The rebuildFrequency should be chosen, s.t. the particles do not move more
   * than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals
   * @param applicableTraversals all applicable traversals
   * @param cellSizeFactor cell size factor relative to cutoff. Verlet lists are only implemented for values >= 1.0
   * (smaller values are set to 1.0).
   */
  VerletListsLinkedBase(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                        const double skin, const unsigned int rebuildFrequency,
                        const std::set<TraversalOption> &applicableTraversals, const double cellSizeFactor)
      : ParticleContainer<Particle, FullParticleCell<Particle>>(boxMin, boxMax, cutoff, skin),
        _linkedCells(boxMin, boxMax, cutoff, skin, std::max(1.0, cellSizeFactor)),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _neighborListIsValid(false),
        _verletBuiltNewton3(false) {
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
  unsigned long getNumParticles() override { return _linkedCells.getNumParticles(); }

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
  std::vector<Particle> AUTOPAS_WARN_UNUSED_RESULT updateContainer() override {
    AutoPasLog(debug, "updating container");
    _neighborListIsValid = false;
    return _linkedCells.updateContainer();
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::isContainerUpdateNeeded()
   */
  bool isContainerUpdateNeeded() override {
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // TODO: find a sensible value for ???
#pragma omp parallel for shared(outlierFound)  // if (this->_cells.size() / omp_get_max_threads() > ???)
#endif
    for (size_t cellIndex1d = 0; cellIndex1d < _linkedCells.getCells().size(); ++cellIndex1d) {
      std::array<double, 3> boxmin{0., 0., 0.};
      std::array<double, 3> boxmax{0., 0., 0.};
      _linkedCells.getCellBlock().getCellBoundingBox(cellIndex1d, boxmin, boxmax);
      boxmin = ArrayMath::subScalar(boxmin, this->getSkin() / 2.);
      boxmax = ArrayMath::addScalar(boxmax, this->getSkin() / 2.);
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
      bool updated = internal::checkParticleInCellAndUpdate(*cellptr, pCopy);
      if (updated) {
        return true;
      }
    }
    AutoPasLog(trace,
               "VerletLists: updateHaloParticle was not able to update particle at "
               "[{}, {}, {}]",
               pCopy.getR()[0], pCopy.getR()[1], pCopy.getR()[2]);
    return false;
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return _linkedCells.begin(behavior);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  ParticleIteratorWrapper<Particle> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                      const std::array<double, 3> &higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    return _linkedCells.getRegionIterator(lowerCorner, higherCorner, behavior, true);
  }

  /**
   * Get the dimension of the used cellblock including the haloboxes.
   * @return the dimensions of the used cellblock
   */
  const std::array<std::size_t, 3> &getCellsPerDimension() {
    return _linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
  }

  /**
   * Specifies whether the neighbor lists need to rebuild when using the given Newton 3 option.
   * @param useNewton3 Specifies if newton3 should be used.
   * @return True if the neighbor lists need to be rebuild, false otherwise.
   */
  bool needsRebuild(bool useNewton3) {
    AutoPasLog(debug, "VerletLists: neighborlist is valid: {}", _neighborListIsValid);
    // if the neighbor list is NOT valid, we have not rebuild for _rebuildFrequency steps or useNewton3 changed
    return (not _neighborListIsValid) or (_traversalsSinceLastRebuild >= _rebuildFrequency) or
           (useNewton3 != _verletBuiltNewton3);
  }

  /**
   * Specifies whether the neighbor lists need to be rebuild.
   * @note Assumes that the newton3 type has NOT changed!
   * @return True if the neighbor lists need to be rebuild, false otherwise.
   */
  bool needsRebuild() { return needsRebuild(_verletBuiltNewton3); }

  /**
   * Generates a traversal selector info for this container.
   * @return Traversal selector info for this container.
   */
  TraversalSelectorInfo<ParticleCell> getTraversalSelectorInfo() override {
    return TraversalSelectorInfo<ParticleCell>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
  }

 protected:
  /// internal linked cells storage, handles Particle storage and used to build verlet lists
  LinkedCells<Particle, LinkedParticleCell, LinkedSoAArraysType> _linkedCells;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  /// specifies if the neighbor list is currently valid
  bool _neighborListIsValid;

  /// specifies if the current verlet list was built for newton3
  bool _verletBuiltNewton3;
};

}  // namespace autopas
