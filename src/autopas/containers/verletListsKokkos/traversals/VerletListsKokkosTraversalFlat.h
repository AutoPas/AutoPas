/**
 * @file VerletListsKokkosTraversalFlat.h
 * @date 28.05.2026
 * @author Franziska Duhr
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/verletListsKokkos/traversals/VerletListsKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

#include <Kokkos_Core.hpp>

namespace autopas {

/**
 * This class defines the traversal typically used by the VerletListsKokkos container
 *
 * @tparam Functor
 * @tparam Particle_T
 */
template <class Functor, class Particle_T>
class VerletListsKokkosTraversalFlat : public TraversalInterface, public VerletListsKokkosTraversalInterface<Particle_T> {

};
} 
