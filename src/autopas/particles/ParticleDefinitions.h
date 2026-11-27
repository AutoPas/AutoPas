/**
 * @file ParticleDefinitions.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "KokkosParticleBase.h"
#include "autopas/particles/ParticleBase.h"

namespace autopas {
/**
 * Particle with all variables in 64 bit precision
 */
using ParticleBaseFP64 = ParticleBase<double, unsigned long>;

/**
 * Particle with all variables in 32 bit precision
 */
using ParticleBaseFP32 = ParticleBase<float, unsigned int>;

#ifdef AUTOPAS_ENABLE_KOKKOS
/**
 * Kokkos-compatible Particle with all variables in 64 bit precision
 */
using KokkosParticleBaseFP64 = KokkosParticleBase<double, unsigned long>;

/**
 * Kokkos-compatible Particle with all variables in 32 bit precision
 */
using KokkosParticleBaseFP32 = KokkosParticleBase<float, unsigned int>;
#endif

}  // namespace autopas
