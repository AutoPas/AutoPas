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
using ParticleFP32 = ParticleBase<float, uint64_t>;
/**
 * Particle with all variables in 64 bit precision
 */
using ParticleFP64 = ParticleBase<double, uint64_t>;
/**
 * Alias for Particle with all variables in 64 bit precision
 */
using Particle = ParticleFP64;

}  // namespace autopas
