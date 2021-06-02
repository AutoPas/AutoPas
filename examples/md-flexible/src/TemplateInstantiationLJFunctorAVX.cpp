
#include "autopas/AutoPasImpl.h"
#include "TypeDefinitions.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"

template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::LJFunctorAVX<ParticleType, true, true>*);