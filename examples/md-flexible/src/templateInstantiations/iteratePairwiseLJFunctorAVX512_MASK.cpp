#if defined(MD_FLEXIBLE_FUNCTOR_AVX512_MASK)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(LJFunctorTypeAVX512_MASK *);
//! @endcond

#endif