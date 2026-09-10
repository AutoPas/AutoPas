/**
 * @file NeighborIdentificationFunctor.h
 * @author S. Newcome
 * @date 22.07.2026
 */

#pragma once

#include <string>

#include "autopas/baseFunctors/InteractionListGeneratorFunctor.h"

namespace autopas {

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
class NeighborIdentificationFunctor : public InteractionListGeneratorFunctor<Particle_T, false> {
 public:
  using typename InteractionListGeneratorFunctor<Particle_T, false>::NeighborListAoSType;

  /**
   * Constructor
   * @param neighborListsAoS Reference to the neighbor list map. When used with AutoPas::computeInteractions, the map
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
  NeighborIdentificationFunctor(NeighborListAoSType &neighborListsAoS, double interactionLength,
                                bool gatherNewton3Lists)
      : InteractionListGeneratorFunctor<Particle_T, false>(neighborListsAoS, interactionLength, gatherNewton3Lists) {}

  std::string getName() override { return "NeighborIdentificationFunctor"; }
};

}  // namespace autopas