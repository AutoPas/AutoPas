/**
 * @file computeInteractionsPythonTriwise.h
 *
 * Declares the explicit template instantiation of computeInteractions() for the adapter of Triwise functors written
 * in Python, so that other translation units do not instantiate it again.
 */
#pragma once

#include "autopas/AutoPasDecl.h"
#include "src/autopasPython/PythonFunctorAdapters.h"

//! @cond Doxygen_Suppress
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(autopasPython::PythonTriwiseAdapter *);
//! @endcond
