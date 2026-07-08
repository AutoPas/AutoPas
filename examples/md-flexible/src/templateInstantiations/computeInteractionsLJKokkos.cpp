/**
 * @file computeInteractionsLJKokkos.cpp
 *
 */

#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

#ifdef AUTOPAS_ENABLE_KOKKOS

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeKokkos *);
//! @endcond

#endif