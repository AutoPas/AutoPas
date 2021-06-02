#include "autopas/molecularDynamics/LJFunctor.h"

#include "TypeDefinitions.h"

#include "autopas/AutoPasImpl.h"

template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>*);
