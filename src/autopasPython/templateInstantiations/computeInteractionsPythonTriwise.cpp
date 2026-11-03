/**
 * @file computeInteractionsPythonTriwise.cpp
 *
 * Explicit template instantiation of computeInteractions() for the adapter of Triwise functors written in Python.
 * The instantiations for the native functors are the ones of md-flexible.
 */
#include "computeInteractionsPythonTriwise.h"

#include "autopas/AutoPasImpl.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<ParticleType>::computeInteractions(autopasPython::PythonTriwiseAdapter *);
//! @endcond
