
#include "autopas/AutoPasImpl.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true> *);
//! @endcond
