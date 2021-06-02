
#include "autopas/AutoPasImpl.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "src/TypeDefinitions.h"

template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::LJFunctorAVX<ParticleType, true, true>*);