/**
 * @file ParticleConceptTest.cpp
 * @author S. Newcome
 * @date 28.07.26
 *
 * Compile time check that MoleculeLJ satisfies the ParticleType concept.
 */

#include "autopas/utils/ParticleConcept.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

static_assert(autopas::utils::ParticleType<mdLib::MoleculeLJ>, "MoleculeLJ does not satisfy the ParticleType concept!");
