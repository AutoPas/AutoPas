
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "TypeDefinitions.h"
#include "autopas/AutoPasImpl.h"

template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::FlopCounterFunctor<ParticleType>*);