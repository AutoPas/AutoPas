/**
 * @file iteratePairwiseLJFunctorGlobals.cpp
 *
 * Contains a explicit template instantiation for the iteratePairwise() method with the LJFunctor w. globals of the main
 * AutoPas class and the particle type used by md-flexible. This is linked into the md-flexible executable to enable the
 * other compilation units to only declare, but not instantiate this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
#include "autopas/molecularDynamics/LJMultisiteFunctor.h"
//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true> *);
//! @endcond
#else
#include "autopas/molecularDynamics/LJFunctor.h"
//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true> *);
//! @endcond
#endif

#endif
