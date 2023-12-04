/**
 * @file computeInteractionsFlopCounterFunctor.cpp
 *
 * Contains a explicit template instantiation for the computeInteractions() method with the FlopCounterFunctor of the
 * main AutoPas class and the particle type used by md-flexible, as determined by whether md-flexible is compiled with
 * Multi-Site support or not. This is linked into the md-flexible executable to enable the other compilation units to
 * only declare, but not instantiate this template.
 */

#include "autopas/AutoPasImpl.h"
#include "autopas/baseFunctors/FlopCounterFunctor.h"
#include "src/TypeDefinitions.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(
    autopas::FlopCounterFunctor<ParticleType, LJFunctorTypeAbstract> *);
//! @endcond
