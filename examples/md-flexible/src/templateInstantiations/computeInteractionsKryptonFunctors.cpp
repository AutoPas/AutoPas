/**
 * @file computeInteractionsKryptonFunctors.cpp
 *
 * Contains explicit template instantiations for the computeInteractions() method with the appropriate ab-initio
 * Krypton pair Functor aswell as the non-additive extended Axilrod-Teller-Muto 3-body functor. This is linked into the
 * md-flexible executable to enable the other compilation units to only declare, but not instantiate this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_KRYPTON)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(KryptonPairwiseFunctorType *);
template bool autopas::AutoPas<ParticleType>::computeInteractions(KryptonTriwiseFunctorType *);
//! @endcond

#endif
