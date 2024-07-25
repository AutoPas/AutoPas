/**
 * @file computeInteractionsATFunctorAVX512.cpp
 *
 * Contains an explicit template instantiation for the computeInteractions() method with the appropriate AVX512 intrinsics
 * Axilrod-Teller Functor and Particle Type, as determined by whether md-flexible is compiled with or without Multi-Site
 * support. This is linked into the md-flexible executable to enable the other compilation units to only declare, but
 * not instantiate this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_AT_AVX512)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(ATFunctorAVX512 *);
//! @endcond

#endif
