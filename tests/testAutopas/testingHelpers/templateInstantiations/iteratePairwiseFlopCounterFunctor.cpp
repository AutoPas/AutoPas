/**
 * @file iteratePairwiseFlopCounterFunctor.cpp
 *
 * Contains a explicit template instantiation for the iteratePairwise() method with the FlopCounterFunctor of the main
 * AutoPas class and the particle type used by md-flexible. This is linked into the md-flexible executable to enable the
 * other compilation units to only declare, but not instantiate this template.
 */

#include "autopas/AutoPasImpl.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "testingHelpers/commonTypedefs.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<Particle>::iteratePairwise(autopas::FlopCounterFunctor<Particle> *);
//! @endcond
