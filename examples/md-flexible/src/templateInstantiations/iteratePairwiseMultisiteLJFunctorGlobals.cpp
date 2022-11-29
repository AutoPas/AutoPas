/**
 * @file iteratePairwiseLJMultisiteFunctor.cpp
 *
 * Contains a explicit template instantiation for the iteratePairwise() method with the LJMultisiteFunctor of the main
 * AutoPas class and the particle type used by md-flexible. This is linked into the md-flexible executable to enable the
 * other compilation units to only declare, but not instantiate this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
#include "autopas/AutoPasImpl.h"
#include "autopas/molecularDynamics/LJMultisiteFunctor.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true, true> *);
//! @endcond
#endif
