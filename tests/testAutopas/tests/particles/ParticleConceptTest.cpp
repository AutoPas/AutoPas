/**
 * @file ParticleConceptTest.cpp
 * @author S. Newcome
 * @date 28.07.26
 *
 * Compile time checks that the ParticleType concept accepts the particles AutoPas ships with and rejects types that
 * do not implement the full interface.
 */

#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/ParticleConcept.h"
#include "testingHelpers/NonConstructibleParticle.h"

static_assert(autopas::utils::ParticleType<autopas::ParticleBaseFP64>,
              "ParticleBaseFP64 does not satisfy the ParticleType concept!");

static_assert(autopas::utils::ParticleType<NonConstructibleParticle>,
              "NonConstructibleParticle does not satisfy the ParticleType concept!");

namespace {
/**
 * A particle that hides one of the required members.
 */
struct ParticleWithoutGetF : public autopas::ParticleBaseFP64 {
 private:
  using autopas::ParticleBaseFP64::getF;
};

/**
 * A type that is not a particle at all.
 */
struct NotAParticle {};
}  // namespace

static_assert(not autopas::utils::ParticleType<ParticleWithoutGetF>,
              "A particle missing getF() must not satisfy the ParticleType concept.");

static_assert(not autopas::utils::ParticleType<NotAParticle>,
              "A type unrelated to particles must not satisfy the ParticleType concept.");
