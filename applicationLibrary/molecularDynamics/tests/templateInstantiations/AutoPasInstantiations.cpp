/**
 * @file AutoPasInstantiations.cpp
 *
 * Contains a explicit template instantiation so that they are not generated on a test-by-test basis.
 * Generating similar templates in the same translation unit also saves memory.
 */

#include "autopas/AutoPasImpl.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

//! @cond Doxygen_Suppress
template class autopas::AutoPas<Molecule>;
// clang-format off
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::LJFunctor<Molecule,
                     /* shifting */ false,
                     /*mixing*/ false,
                     autopas::FunctorN3Modes::Both,
                     /*globals*/ false,
                     /*countFLOPs*/ false,
                     /*relevantForTuning*/ true> *);
// clang-format on
//! @endcond
