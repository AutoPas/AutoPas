/**
 * @file Particle.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/particles/ParticleBase.h"

namespace autopas {

/**
 * Particle with all variables in 32 bit precision
 */
using ParticleFP32 = ParticleBase<float, float, unsigned long>;
/**
 * Particle with all variables in 64 bit precision
 */
using ParticleFP64 = ParticleBase<double, double, unsigned long>;
/**
 * Particle with mixed precision for calculations and accumulation
 */
using ParticleMP = ParticleBase<float, double, unsigned long>;

#if AUTOPAS_PRECISION_MODE == SPSP
using Particle = ParticleFP32;
#elif AUTOPAS_PRECISION_MODE == DPDP
using Particle = ParticleFP64;
#else
using Particle = ParticleMP;
#endif

}  // namespace autopas
