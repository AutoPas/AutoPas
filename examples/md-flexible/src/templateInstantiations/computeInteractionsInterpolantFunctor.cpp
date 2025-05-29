/**
 * @file computeInteractionsLJInterpolantFunctor.cpp
 *
 * Contains an explicit template instantiation for the computeInteractions() method with the appropriate Interpolated
 * Functors and Particle Type, as determined by whether md-flexible is compiled with or without Multi-Site support. This
 * is linked into the md-flexible executable to enable the other compilation units to only declare, but not instantiate
 * this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_PAIRWISE_INTERPOLANT)

#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(LJInterpolantFunctorType *);
//! @endcond

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(KryptonPairInterpolantFunctorType *);
//! @endcond

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(ArgonPairInterpolantFunctorType *);
//! @endcond
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_TRIWISE_INTERPOLANT)

#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(KryptonTripletInterpolantFunctorType *);
//! @endcond

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(ArgonTripletInterpolantFunctorType *);
//! @endcond
#endif