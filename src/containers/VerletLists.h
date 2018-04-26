/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "LinkedCells.h"
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
  // AOS
  typedef std::map<Particle*, std::vector<Particle*>>
      AoS_verletlist_storage_type;

  // SOA
  typedef std::map<decltype(Particle().getID()), size_t>
      particleid_to_verletlistindex_container_type;
  typedef std::vector<decltype(Particle().getID())>
      verletlistindex_to_particleid_container_type;

 public:
  /**
   * Constructor of the VerletLists class.
   * Each cell of the verlet lists class, is at least of size cutoff + skin.
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
    if (needsRebuild()) {  // if we need to rebuild the list, we should rebuild it!
      this->updateVerletListsAoS(useNewton3);
    }
    this->iterateVerletListsAoS(f, useNewton3);
    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  AoS_verletlist_storage_type& getVerletListsAoS() { return _verletListsAoS; }

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
  void updateContainer() override{
    _neighborListIsValid = false;
    LinkedCells<Particle, ParticleCell>::updateContainer();
  }

 private:
  class VerletListGeneratorFunctor
      : public autopas::Functor<Particle, ParticleCell> {
   public:
    VerletListGeneratorFunctor(AoS_verletlist_storage_type& verletListsAoS,
                               particleid_to_verletlistindex_container_type&
                                   particleIDtoVerletListIndexMap,
                               double cutoffskinsquared)
        : _verletListsAoS(verletListsAoS),
          _particleIDtoVerletListIndexMap(particleIDtoVerletListIndexMap),
          _cutoffskinsquared(cutoffskinsquared) {}

    void AoSFunctor(Particle& i, Particle& j, bool newton3 = true) override {
      auto dist = arrayMath::sub(i.getR(), j.getR());
      double distsquare = arrayMath::dot(dist, dist);
      if (distsquare < _cutoffskinsquared)
        // this is thread safe, only if particle i is accessed by only one
        // thread at a time. which is ensured, as particle i resides in a
        // specific cell and each cell is only accessed by one thread at a time
        // (ensured by traversals)
        // also the list is not allowed to be resized!

        _verletListsAoS[&i].push_back(&j);
    }

   private:
    AoS_verletlist_storage_type& _verletListsAoS;
    particleid_to_verletlistindex_container_type&
        _particleIDtoVerletListIndexMap;
    double _cutoffskinsquared;
  };

  bool needsRebuild() {
    return (not _neighborListIsValid)  // if the neighborlist is NOT valid a
                                       // rebuild is needed
           or (_traversalsSinceLastRebuild >=
               _rebuildFrequency);  // rebuild with frequency
  }

  void updateVerletListsAoS(bool useNewton3) {
    _verletListsAoS.clear();
    updateIdMapAoS();
    VerletListGeneratorFunctor f(_verletListsAoS,
                                 _particleIDtoVerletListIndexContainer,
                                 (this->getCutoff() * this->getCutoff()));

    LinkedCells<Particle, ParticleCell>::iteratePairwiseAoS2(&f, useNewton3);

    // the neighbor list is now valid
    _neighborListIsValid = true;
  }

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

  size_t updateIdMapAoS() {
    size_t i = 0;

    // DON'T simply parallelize this loop!!! this needs modifications if you
    // want to parallelize it!
    for (auto iter = this->begin(); iter.isValid(); ++iter, ++i) {
      _verletListsAoS[&(*iter)];
    }

    return i;
  }

  /// map that converts the id type of the particle to the actual position in
  /// the verlet list.
  /// This is needed, as the particles don't have a local id.
  /// @todo remove this and add local id to the particles (this saves a
  /// relatively costly lookup)
  particleid_to_verletlistindex_container_type
      _particleIDtoVerletListIndexContainer;

  verletlistindex_to_particleid_container_type
      _verletListIndextoParticleIDContainer;

  /// verlet lists.
  AoS_verletlist_storage_type _verletListsAoS;

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
