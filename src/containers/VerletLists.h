/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "LinkedCells.h"
#include "VerletListHelpers.h"
#include "utils/arrayMath.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise
 * interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @note This class does NOT work with RMM cells and is not intended to!
 * @tparam Particle
 * @tparam ParticleCell
 * @todo deleting particles should also invalidate the verlet lists - should be
 * implemented somehow
 */
template <class Particle, class ParticleCell>
class VerletLists : public LinkedCells<Particle, ParticleCell> {
  typedef VerletListHelpers<Particle, ParticleCell> verlet_internal;

 public:
  /**
   * Constructor of the VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * The rebuildFrequency should be chosen, s.t. the particles do not move more
   * than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals.
   */
  VerletLists(const std::array<double, 3> boxMin,
              const std::array<double, 3> boxMax, double cutoff, double skin,
              unsigned int rebuildFrequency = 1)
      : LinkedCells<Particle, ParticleCell>(boxMin, boxMax, cutoff + skin),
        _skin(skin),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _neighborListIsValid(false) {}

  void iteratePairwiseAoS(Functor<Particle, ParticleCell>* f,
                          bool useNewton3 = true) override {
    iteratePairwiseAoS2(f, useNewton3);
  }

  /**
   * same as iteratePairwiseAoS, but potentially faster (if called with the
   * derived functor), as the class of the functor is known and thus the
   * compiler can do some better optimizations.
   *
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3 defines whether newton3 should be used
   */
  template <class ParticleFunctor>
  void iteratePairwiseAoS2(ParticleFunctor* f, bool useNewton3 = true) {
    if (needsRebuild()) {  // if we need to rebuild the list, we should rebuild
                           // it!
      this->updateVerletListsAoS(useNewton3);
      // the neighbor list is now valid
      _neighborListIsValid = true;
      _traversalsSinceLastRebuild = 0;
    }
    this->iterateVerletListsAoS(f, useNewton3);
    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  typename verlet_internal::AoS_verletlist_storage_type& getVerletListsAoS() {
    return _verletListsAoS;
  }

  /**
   * @copydoc LinkedCells::addParticle()
   * @note this function invalidates the neighbor lists
   */
  void addParticle(Particle& p) override {
    _neighborListIsValid = false;
    LinkedCells<Particle, ParticleCell>::addParticle(p);
  }

  /**
   * @copydoc LinkedCells::addHaloParticle()
   * @note this function invalidates the neighbor lists
   */
  void addHaloParticle(Particle& haloParticle) override {
    _neighborListIsValid = false;
    LinkedCells<Particle, ParticleCell>::addHaloParticle(haloParticle);
  }

  /**
   * @copydoc LinkedCells::updateContainer()
   * @note this function invalidates the neighbor lists
   */
  void updateContainer() override {
    _neighborListIsValid = false;
    LinkedCells<Particle, ParticleCell>::updateContainer();
  }

  /**
   * Checks whether the neighbor lists are valid.
   * A neighbor list is valid if all pairs of particles whose interaction should
   * be calculated are represented in the neighbor lists.
   * @return whether the list is valid
   * @note this check involves pair-wise interaction checks and is thus
   * relatively costly.
   */
  bool checkNeighborListsAreValid(bool useNewton3 = true) {
    // if a particle was added or deleted, ... the list is definitely invalid
    if (not _neighborListIsValid) {
      return false;
    }
    // if a particle moved more than skin/2 outside of its cell the list is
    // invalid
    if (this->isContainerUpdateNeeded()) {
      return false;
    }

    // particles can also simply be very close already:
    typename verlet_internal::VerletListValidityCheckerFunctor
        validityCheckerFunctor(
            _verletListsAoS, _particleIDtoVerletListIndexContainer,
            ((this->getCutoff() - _skin) * (this->getCutoff() - _skin)));

    LinkedCells<Particle, ParticleCell>::iteratePairwiseAoS2(
        &validityCheckerFunctor, useNewton3);

    return validityCheckerFunctor.neighborlistsAreValid();
  }

  bool isContainerUpdateNeeded() override {
    for (int cellIndex1d = 0; cellIndex1d < this->_data.size(); ++cellIndex1d) {
      std::array<double, 3> boxmin;
      std::array<double, 3> boxmax;
      this->_cellBlock.getCellBoundingBox(cellIndex1d, boxmin, boxmax);
      boxmin = arrayMath::addScalar(boxmin, -_skin / 2.);
      boxmax = arrayMath::addScalar(boxmax, +_skin / 2.);
      for (auto iter = this->_data[cellIndex1d].begin(); iter.isValid();
           ++iter) {
        if (not iter->inBox(boxmin, boxmax)) {
          return true;  // we need an update
        }
      }
    }
    return false;
  }

 protected:
  /**
   * specifies whether the neighbor lists need to be rebuild
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild() {
    AutoPasLogger->debug("VerletLists: neighborlist is valid: {}",_neighborListIsValid);
    return (not _neighborListIsValid)  // if the neighborlist is NOT valid a
                                       // rebuild is needed
           or (_traversalsSinceLastRebuild >=
               _rebuildFrequency);  // rebuild with frequency
  }

  /**
   * update the verlet lists for AoS usage
   * @param useNewton3
   */
  virtual void updateVerletListsAoS(bool useNewton3) {
    _verletListsAoS.clear();
    updateIdMapAoS();
    typename verlet_internal::VerletListGeneratorFunctor f(
        _verletListsAoS, _particleIDtoVerletListIndexContainer,
        (this->getCutoff() * this->getCutoff()));

    LinkedCells<Particle, ParticleCell>::iteratePairwiseAoS2(&f, useNewton3);
  }

  /**
   * iterate over the verlet lists using the AoS traversal
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iterateVerletListsAoS(ParticleFunctor* f, const bool useNewton3) {
    /// @todo optimize iterateVerletListsAoS, e.g. by using traversals with
    /// different

    // don't parallelize this with a simple openmp, unless useNewton3=false
    for (auto& list : _verletListsAoS) {
      Particle& i = *list.first;
      for (auto j_ptr : list.second) {
        Particle& j = *j_ptr;
        f->AoSFunctor(i, j, useNewton3);
      }
    }
  }

  /**
   * update the AoS id maps.
   * The Id Map is used to map the id of a particle to the actual particle
   * @return
   */
  size_t updateIdMapAoS() {
    size_t i = 0;

    // DON'T simply parallelize this loop!!! this needs modifications if you
    // want to parallelize it!
    for (auto iter = this->begin(); iter.isValid(); ++iter, ++i) {
      _verletListsAoS[&(*iter)];
    }

    return i;
  }

 private:
  /// map that converts the id type of the particle to the actual position in
  /// the verlet list.
  /// This is needed, as the particles don't have a local id.
  /// @todo remove this and add local id to the particles (this saves a
  /// relatively costly lookup)
  typename verlet_internal::particleid_to_verletlistindex_container_type
      _particleIDtoVerletListIndexContainer;

  typename verlet_internal::verletlistindex_to_particleid_container_type
      _verletListIndextoParticleIDContainer;

  /// verlet lists.
  typename verlet_internal::AoS_verletlist_storage_type _verletListsAoS;

  double _skin;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  // specifies if the neighbor list is currently valid
  bool _neighborListIsValid;
};

} /* namespace autopas */
