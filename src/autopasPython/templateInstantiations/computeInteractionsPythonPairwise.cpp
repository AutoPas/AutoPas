/**
 * @file computeInteractionsPythonPairwise.cpp
 *
 * Explicit template instantiation of computeInteractions() for the adapter of Pairwise functors written in Python.
 * The instantiations for the native functors are the ones of md-flexible.
 */
#include "computeInteractionsPythonPairwise.h"

#include "autopas/AutoPasImpl.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(autopasPython::PythonPairwiseAdapter *);
//! @endcond
