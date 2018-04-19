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
 * This class builds neighbour lists for the particle interactions.
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class VerletLists : public LinkedCells<Particle, ParticleCell> {
  typedef std::vector<std::vector<size_t>> verletlist_storage_type;
  typedef std::map<decltype(Particle().getID()), size_t>
      particleid_to_verletlistindex_container_type;
  typedef std::vector<decltype(Particle().getID())>
      verletlistindex_to_particleid_container_type;

 public:
  /**
   * Constructor of the VerletLists class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
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
    this->updateVerletLists(useNewton3);
    this->iterateVerletListsAoS(f, useNewton3);
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  verletlist_storage_type& getVerletLists() { return _verletLists; }

 private:
  class VerletListGeneratorFunctor
      : public autopas::Functor<Particle, ParticleCell> {
   public:
    VerletListGeneratorFunctor(
        verletlist_storage_type& verletLists,
        particleid_to_verletlistindex_container_type& particleIDtoVerletListIndexMap,
        double cutoffskinsquared)
        : _verletLists(verletLists),
          _particleIDtoVerletListIndexMap(particleIDtoVerletListIndexMap),
          _cutoffskinsquared(cutoffskinsquared) {}

    void AoSFunctor(Particle& i, Particle& j, bool newton3 = true) override {
      auto dist = arrayMath::sub(i.getR(), j.getR());
      double distsquare = arrayMath::dot(dist, dist);

      if (distsquare < _cutoffskinsquared)
        _verletLists[_particleIDtoVerletListIndexMap[i.getID()]].push_back(
            _particleIDtoVerletListIndexMap[j.getID()]);
    }

   private:
    verletlist_storage_type& _verletLists;
    particleid_to_verletlistindex_container_type& _particleIDtoVerletListIndexMap;
    double _cutoffskinsquared;
  };

  void updateVerletLists(bool useNewton3) {
    size_t particleNumber = updateIdMap();
    _verletLists.clear();
    _verletLists.resize(particleNumber);
    VerletListGeneratorFunctor f(_verletLists, _particleIDtoVerletListIndexContainer,
                                 (this->getCutoff() * this->getCutoff()));

    LinkedCells<Particle, ParticleCell>::iteratePairwiseAoS2(&f, useNewton3);
  }

  template <class ParticleFunctor>
  void iterateVerletListsAoS(ParticleFunctor* f, const bool useNewton3) {
    for(auto& list : _verletLists){

    }
  }

  size_t updateIdMap() {
    _particleIDtoVerletListIndexContainer.clear();  // this is not necessary, but it
                                              // does not hurt too much,
                                              // probably...
    size_t i = 0;
    for (auto iter = this->begin(); iter.isValid(); ++iter, ++i) {
      _particleIDtoVerletListIndexContainer[iter->getID()] = i;
      _verletListIndextoParticleIDContainer.push_back(iter->getID());
    }

    return i;
  }

  /// map that converts the id type of the particle to the actual position in
  /// the verlet list.
  /// This is needed, as the particles don't have a local id.
  /// @todo remove this and add local id to the particles (this saves a
  /// relatively costly lookup)
  particleid_to_verletlistindex_container_type _particleIDtoVerletListIndexContainer;

  verletlistindex_to_particleid_container_type _verletListIndextoParticleIDContainer;

  /// verlet lists.
  verletlist_storage_type _verletLists;

  double _skin;
};

} /* namespace autopas */
