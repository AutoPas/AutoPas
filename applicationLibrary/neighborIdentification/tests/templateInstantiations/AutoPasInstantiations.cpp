/**
 * @file AutoPasInstantiations.cpp
 *
 * Contains a explicit template instantiation so that they are not generated on a test-by-test basis.
 * Generating similar templates in the same translation unit also saves memory.
 */

#include "autopas/AutoPasImpl.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "neighborIdentificationLibrary/NeighborIdentificationFunctor.h"

//! @cond Doxygen_Suppress
template class autopas::AutoPas<autopas::ParticleBaseFP64>;
template bool autopas::AutoPas<autopas::ParticleBaseFP64>::computeInteractions(
    autopas::NeighborIdentificationFunctor<autopas::ParticleBaseFP64> *);
//! @endcond
