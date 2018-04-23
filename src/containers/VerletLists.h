/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "LinkedCells.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @todo This class does not yet provide the ability to reuse computet verlet lists.
 * @note This class does NOT work with RMM cells and is not intended to!
 * @tparam Particle
 * @tparam ParticleCell
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
   *
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   */
  VerletLists(const std::array<double, 3> boxMin,
              const std::array<double, 3> boxMax, double cutoff, double skin)
      : LinkedCells<Particle, ParticleCell>(boxMin, boxMax, cutoff + skin),
        _skin(skin) {}

  void iteratePairwiseAoS(Functor<Particle, ParticleCell>* f,
                          bool useNewton3 = true) override {
    iteratePairwiseAoS2(f, useNewton3);
  }

  /**
   * same as iteratePairwiseAoS, but potentially faster (if called with the
   * derived functor), as the class of the functor is known and thus the
   * compiler can do some better optimizations.
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3 defines whether newton3 should be used
   */
  template <class ParticleFunctor>
  void iteratePairwiseAoS2(ParticleFunctor* f, bool useNewton3 = true) {
    this->updateVerletListsAoS(useNewton3);
    this->iterateVerletListsAoS(f, useNewton3);
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  AoS_verletlist_storage_type& getVerletListsAoS() { return _verletListsAoS; }

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

  void updateVerletListsAoS(bool useNewton3) {
    _verletListsAoS.clear();
    updateIdMapAoS();
    VerletListGeneratorFunctor f(_verletListsAoS,
                                 _particleIDtoVerletListIndexContainer,
                                 (this->getCutoff() * this->getCutoff()));

    LinkedCells<Particle, ParticleCell>::iteratePairwiseAoS2(&f, useNewton3);
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
};

} /* namespace autopas */
