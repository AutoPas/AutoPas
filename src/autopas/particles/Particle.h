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
using ParticleFP32 = ParticleBase<float, unsigned long>;
/**
 * Particle with all variables in 64 bit precision
 */
using ParticleFP64 = ParticleBase<double, unsigned long>;
/**
 * Alias for Particle with unsigned long id
 */
template<typename floatType>
using Particle = ParticleBase<floatType, unsigned long>;

}  // namespace autopas
