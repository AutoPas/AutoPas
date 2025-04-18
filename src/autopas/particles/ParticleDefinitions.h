/**
 * @file ParticleDefinitions.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/particles/ParticleBase.h"

namespace autopas {
/**
 * Particle with all variables in 64 bit precision
 */
using ParticleBaseFP64 = ParticleBase<double, unsigned long>;

}  // namespace autopas
