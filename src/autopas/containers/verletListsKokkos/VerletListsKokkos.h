/**
 * @file VerletListsKokkos.h
 * @date 28.05.2026
 * @author Franziska Duhr
 */
#pragma once

#include "autopas/containers/ParticleContainerInterface.h"

namespace autopas {

/**
 * @class VerletListsKokkos
 * @brief A container for managing particles using Kokkos for parallel computation.
 * @tparam Particle_T AoS Particle type that is used with this container
 */
template <class Particle_T>
class VerletListsKokkos : public ParticleContainerInterface<Particle_T> {
};

} 
