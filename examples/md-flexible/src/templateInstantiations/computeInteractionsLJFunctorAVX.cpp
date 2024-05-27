/**
 * @file computeInteractionsLJFunctorAVX.cpp
 *
 * Contains an explicit template instantiation for the computeInteractions() method with the appropriate AVX-vectorized
 * Lennard-Jones Functor and Particle Type, as determined by whether md-flexible is compiled with or without Multi-Site
 * support. This is linked into the md-flexible executable to enable the other compilation units to only declare, but
 * not instantiate this template.
 */
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeAVX *);
//! @endcond

#endif