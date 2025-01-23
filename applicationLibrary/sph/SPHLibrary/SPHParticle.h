/**
 * @file SPHParticle.h
 *
 * @date 21/Jan/2025
 * @author huberf
 */

#pragma once

#include "SPHParticleBase.h"

namespace sphLib {
/**
 * SPHParticle with all variables in 32 bit precision
 */
using SPHParticleFP32 = SPHParticleBase<float, float, unsigned long>;
/**
 * SPHParticle with all variables in 64 bit precision
 */
using SPHParticleFP64 = SPHParticleBase<double, double, unsigned long>;
/**
 * SPHParticle with mixed precision for calculations and accumulation
 */
using SPHParticleMP = SPHParticleBase<float, double, unsigned long>;

#if AUTOPAS_PRECISION_MODE == SPSP
using SPHParticle = SPHParticleFP32;
#elif AUTOPAS_PRECISION_MODE == SPDP
using SPHParticle = SPHParticleMP;
#else
using SPHParticle = SPHParticleFP64;
#endif
}  // namespace sphLib