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
typedef ParticleBase<float, unsigned long> ParticleFP32;
/**
 * Particle with all variables in 64 bit precision
 */
typedef ParticleBase<double, unsigned long> ParticleFP64;
/**
 * Alias for Particle with all variables in 64 bit precision
 */
typedef ParticleFP64 Particle;


}  // namespace autopas
