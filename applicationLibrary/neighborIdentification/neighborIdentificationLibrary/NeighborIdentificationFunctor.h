/**
 * @file NeighborIdentificationFunctor.h
 * @author S. Newcome
 * @date 22.07.2026
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "autopas/baseFunctors/InteractionListGeneratorFunctor.h"

namespace autopas {

/**
 * A policy for managing neighbor lists in an Array of Structures (AoS) layout.
 * This is needed for the NeighborIdentificationFunctor.
 * @tparam Particle_T The type of Particle class used.
 */
template <class Particle_T>
class AoSNeighborListPolicy {
 public:
  /**
   * Constructor.
   * @param neighborListAoS neighbor list map to be filled by the functor.
   */
  explicit AoSNeighborListPolicy(std::unordered_map<Particle_T *, std::vector<Particle_T *>> &neighborListAoS)
      : _neighborListsAoS(neighborListAoS) {}

  /**
   * Initializes the neighbor list map for the given range of particles: clears any previous contents and inserts an
   * empty neighbor list for every particle. This must happen before the functor is applied, and it must also happen
   * after any particle rearrangement in memory.
   *
   * @tparam ParticleIterator_T Type of the particle iterator
   * @param particlesBegin Iterator to the first particle. Iterated until no longer valid.
   */
  template <class ParticleIterator_T>
  void initializeNeighborList(ParticleIterator_T particlesBegin) {
    _neighborListsAoS.clear();
    for (auto iter = particlesBegin; iter.isValid(); ++iter) {
      _neighborListsAoS[&(*iter)];
    }
  }

  /**
   * Adds a neighbor j to the neighbor list of particle i.
   * @param i The first particle.
   * @param j The neighbor particle.
   */
  void add(Particle_T *i, Particle_T *j) { _neighborListsAoS[i].push_back(j); }

 private:
  std::unordered_map<Particle_T *, std::vector<Particle_T *>> &_neighborListsAoS;
};

/**
 * This functor generates lists of particles within interactionLength of each other, which could be used within
 * user simulators to replace their neighbor identification/contact detection/cutoff check passes, although they should
 * consider the warning below.
 *
 * In some fields this can be called "cutoff checking" or "broad-phase contact detection" but we will refer to this
 * in this file as neighbor identification.
 *
 * @warning To take fully advantage of AutoPas, we recommend, instead of using this functor and apply forces externally
 * to AutoPas, using a functor which  integrates neighbor identification and force calculation in one call over using
 * this functor (see @ref LJFunctor.h as an example). We provide this functor primarily for codes which already
 * perform separate neighbor identification and force calculation passes, to easily experience some benefit of AutoPas
 * with minimal code refactors.
 *
 * @details After applying AutoPas's computeInteractions function with this functor, a
 * std::unordered_map<Particle_T *, std::vector<Particle_T *>> is filled, mapping from each particle pointer to a vector
 * of pointers to all particles within interactionLength of the first.
 *
 * @tparam Particle_T The type of Particle class used.
 */
template <class Particle_T>
class NeighborIdentificationFunctor
    : public InteractionListGeneratorFunctor<Particle_T, AoSNeighborListPolicy<Particle_T>> {
 public:
  /**
   * Constructor
   * @param policy Reference to the neighbor list map. When used with AutoPas::computeInteractions, the map
   * will be overridden.
   * @param interactionLength The distance between particles within which particle pairs get added to the neighbor
   * lists.
   * @param gatherNewton3Lists If false, for a particle pair i, j that are neighbors, **both** particle j will be in
   * particle i's list **and** particle i will be in particle j's list. If true, for a particle pair i, j that are
   * neighbors, **either** particle j will be in particle i's list **or** particle i will be in particle j's list. The
   * latter can be used more easily to reduce calculations by applying forces to both particles in the pair, but care
   * should be taken to avoid race conditions. If true, we make no guarantee which of particle i or j will be in the
   * other's list. If true, this functor will only allow functor calls with newton3 enabled.
   */
  NeighborIdentificationFunctor(AoSNeighborListPolicy<Particle_T> &policy, double interactionLength,
                                const bool gatherNewton3Lists)
      : InteractionListGeneratorFunctor<Particle_T, AoSNeighborListPolicy<Particle_T>>(policy, interactionLength,
                                                                                       gatherNewton3Lists),
        _policy(policy) {}

  std::string getName() override { return "NeighborIdentificationFunctor"; }

  /**
   * Initializes the neighbor list map for the given range of particles: clears any previous contents and inserts an
   * empty neighbor list for every particle. This must happen before the functor is applied, and it must also happen
   * after any particle rearrangement in memory. AutoPas's computeInteractions() calls
   * this automatically when this functor (or a child of it) is used, so the handling of the list is safe when used
   * with computeInteractions.
   *
   * @tparam ParticleIterator_T Type of the particle iterator
   * @param particlesBegin Iterator to the first particle. Iterated until no longer valid.
   */
  template <class ParticleIterator_T>
  void initializeNeighborList(ParticleIterator_T particlesBegin) {
    _policy.initializeNeighborList(particlesBegin);
  }

  /**
   * Whether the functor is relevant for tuning.
   * @return true
   */
  bool isRelevantForTuning() override { return true; }

 private:
  AoSNeighborListPolicy<Particle_T> &_policy;
};

}  // namespace autopas