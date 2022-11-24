/**
 * @file iteratePairwiseLJFunctorSVE.cpp
 *
 * Contains a explicit template instantiation for the iteratePairwise() method with the LJFunctorSVE of the main
 * AutoPas class and the particle type used by md-flexible. This is linked into the md-flexible executable to enable the
 * other compilation units to only declare, but not instantiate this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
#include "autopas/AutoPasImpl.h"
#include "autopas/molecularDynamics/LJFunctorSVE.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::LJFunctorSVE<ParticleType, true, true> *);
//! @endcond
#endif
