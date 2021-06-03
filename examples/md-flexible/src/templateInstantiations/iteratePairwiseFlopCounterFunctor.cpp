
#include "autopas/AutoPasImpl.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::FlopCounterFunctor<ParticleType> *);
//! @endcond
