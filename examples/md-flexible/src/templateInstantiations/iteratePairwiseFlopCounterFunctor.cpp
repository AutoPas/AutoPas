/**
 * @file iteratePairwiseFlopCounterFunctor.cpp
 *
 * Contains a explicit template instantiation for the iteratePairwise() method with the FlopCounterFunctor of the main
 * AutoPas class and the particle type used by md-flexible, as determined by whether AutoPas is compiled with Multi-Site support or not.
 * This is linked into the md-flexible executable to enable the other compilation units to only declare, but not instantiate this template.
 */

#include "autopas/AutoPasImpl.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::FlopCounterFunctor<ParticleType> *);
//! @endcond
