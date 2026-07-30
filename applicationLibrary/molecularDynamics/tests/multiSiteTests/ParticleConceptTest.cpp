/**
 * @file ParticleConceptTest.cpp
 * @author S. Newcome
 * @date 28.07.26
 *
 * Compile time check that MultisiteMoleculeLJ satisfies the ParticleType concept.
 */

#include "autopas/utils/ParticleConcept.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"

static_assert(autopas::utils::ParticleType<mdLib::MultisiteMoleculeLJ>,
              "MultisiteMoleculeLJ does not satisfy the ParticleType concept!");
