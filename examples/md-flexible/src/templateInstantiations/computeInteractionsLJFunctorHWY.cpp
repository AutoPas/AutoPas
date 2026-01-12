/**
 * @file computeInteractionsLJFunctorHWY.cpp
 *
 */
#if defined(MD_FLEXIBLE_FUNCTOR_HWY)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeHWY *);
//! @endcond

#endif