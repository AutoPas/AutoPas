/**
 * @file ParticleConceptTest.cpp
 * @author S. Newcome
 * @date 28.07.26
 *
 * Compile time check that SPHParticle satisfies the ParticleType concept.
 */

#include "SPHLibrary/SPHParticle.h"
#include "autopas/utils/ParticleConcept.h"

static_assert(autopas::utils::ParticleType<sphLib::SPHParticle>,
              "SPHParticle does not satisfy the ParticleType concept!");
