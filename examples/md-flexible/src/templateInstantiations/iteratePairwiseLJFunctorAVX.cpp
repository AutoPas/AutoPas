
#include "autopas/AutoPasImpl.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::LJFunctorAVX<ParticleType, true, true> *);
//! @endcond
