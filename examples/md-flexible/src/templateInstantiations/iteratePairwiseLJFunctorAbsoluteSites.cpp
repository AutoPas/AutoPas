//
// Created by johnny on 14.10.23.
//

/**
 * @file iteratePairwiseLJFunctorAbsoluteSites.cpp
 *
 */

#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(LJFunctorTypeAbsPos *);
//! @endcond

#endif
