/**
 * @file VerletListHelpers.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>
namespace autopas {
template <class Particle, class ParticleCell>
/**
 * class of helpers for verlet lists
 * @tparam Particle
 * @tparam ParticleCell
 */
class VerletListHelpers {
 public:
  /// AOS verlet list storage
  typedef std::map<Particle *, std::vector<Particle *>>
      AoS_verletlist_storage_type;

  /**
   * This functor can generate verlet lists using the typical pairwise
   * traversal.
   * @todo: SoA?
   */
  class VerletListGeneratorFunctor
      : public autopas::Functor<Particle, ParticleCell> {
   public:
    /**
     * Constructor
     * @param verletListsAoS
     * @param particleIDtoVerletListIndexMap
     * @param cutoffskinsquared
     */
    VerletListGeneratorFunctor(AoS_verletlist_storage_type &verletListsAoS,
                               double cutoffskinsquared)
        : _verletListsAoS(verletListsAoS),
          _cutoffskinsquared(cutoffskinsquared) {}

    void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
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
    AoS_verletlist_storage_type &_verletListsAoS;
    double _cutoffskinsquared;
  };

  /**
   * This functor checks the validity of neighborhood lists.
   * If a pair of particles has a distance of less than the cutoff radius it
   * checks whether the pair is represented in the verlet list.
   * If the pair is not present in the list the neigborhood lists are invalid
   * and neighborlistsAreValid()  will return false.
   * @todo: SoA?
   */
  class VerletListValidityCheckerFunctor
      : public autopas::Functor<Particle, ParticleCell> {
   public:
    /**
     * Constructor
     * @param verletListsAoS
     * @param particleIDtoVerletListIndexMap
     * @param cutoffsquared
     */
    VerletListValidityCheckerFunctor(
        AoS_verletlist_storage_type &verletListsAoS,
        double cutoffsquared)
        : _verletListsAoS(verletListsAoS),
          _cutoffsquared(cutoffsquared),
          _valid(true) {}

    void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
      auto dist = arrayMath::sub(i.getR(), j.getR());
      double distsquare = arrayMath::dot(dist, dist);
      if (distsquare < _cutoffsquared) {
        // this is thread safe, we have variables on the stack
        auto found = std::find(_verletListsAoS[&i].begin(),
                               _verletListsAoS[&i].end(), &j);
        if (found == _verletListsAoS[&i].end()) {
          // this is thread safe, as _valid is atomic
          _valid = false;
        }
      }
    }

    /**
     * Returns whether the neighbour list are valid.
     * Call this after performing the pairwise traversal
     * @return
     */
    bool neighborlistsAreValid() { return _valid; }

   private:
    AoS_verletlist_storage_type &_verletListsAoS;
    double _cutoffsquared;

    // needs to be thread safe
    std::atomic<bool> _valid;
  };

};  // class verlet_internal
}  // namespace autopas