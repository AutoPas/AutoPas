
#include "autopas/AutoPasImpl.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "src/TypeDefinitions.h"

template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::FlopCounterFunctor<ParticleType>*);