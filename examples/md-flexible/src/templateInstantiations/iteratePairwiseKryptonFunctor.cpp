/**
 * @file iteratePairwiseKryptonFunctor.cpp
 *
 * Contains an explicit template instantiation for the iteratePairwise() method with the appropriate auto-vectorized
 * ab-initio Krypton Functor. This is linked into the md-flexible executable to enable the other compilation units to only declare, but
 * not instantiate this template.
 */

#if defined(MD_FLEXIBLE_FUNCTOR_KRYPTON)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(KryptonFunctorType *);
//! @endcond

#endif
