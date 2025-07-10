/**
 * @file computeInteractionsKryptonFunctor.cpp
 *
 */

#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(KryptonPairAbInitioFunctorType *);
//! @endcond

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(KryptonExtendedATMFunctorType *);
//! @endcond