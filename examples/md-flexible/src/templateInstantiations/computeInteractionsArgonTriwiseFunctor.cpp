/**
* @file computeInteractionsArgonPairwiseFunctor.cpp
*
* Contains an explicit template instantiation for the computeInteractions() method with the appropriate ab-initio
* Argon Triwise Functor. This is linked into the md-flexible executable to enable the other compilation units to only
* declare, but not instantiate this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_ARGON_TRIWISE)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(ArgonTriwiseFunctorType *);
//! @endcond

#endif
