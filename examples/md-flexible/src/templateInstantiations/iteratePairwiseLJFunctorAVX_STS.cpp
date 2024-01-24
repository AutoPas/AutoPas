/**
 * @file iteratePairwiseLJFunctorAVX_STS.cpp
 *
 * Contains an explicit template instantiation for the iteratePairwise() method with the appropriate AVX-vectorized
 * Lennard-Jones Functor (using Site-to-Site cutoffs) and Particle Type. This is linked into the md-flexible executable
 * to enable the other compilation units to only declare, but
 * not instantiate this template.
 */
#if defined(MD_FLEXIBLE_FUNCTOR_AVX_STS) && defined(__AVX__)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(LJFunctorTypeAVXSTS *);
//! @endcond

#endif