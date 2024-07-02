/**
 * @file AutoPasInstantiations.cpp
 *
 * Contains a explicit template instantiation so that they are not generated on a test-by-test basis.
 * Generating similar templates in the same translation unit also saves memory.
 */

#include "autopas/AutoPasImpl.h"
#include "autopas/baseFunctors/FlopCounterFunctor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/NonConstructibleParticle.h"
#include "testingHelpers/commonTypedefs.h"

//! @cond Doxygen_Suppress
template class autopas::AutoPas<Molecule>;
template class autopas::AutoPas<NonConstructibleParticle>;

template bool autopas::AutoPas<Molecule>::computeInteractions(mdLib::LJFunctor<Molecule> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::LJFunctor<Molecule, /* shifting */ true, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                     /*globals*/ true> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::LJFunctor<Molecule, /* shifting */ true, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                     /*globals*/ false> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(EmptyPairwiseFunctor<Molecule> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(
    autopas::FlopCounterFunctor<Molecule, mdLib::LJFunctor<Molecule>> *);
template bool autopas::AutoPas<NonConstructibleParticle>::computeInteractions(
    MockPairwiseFunctor<NonConstructibleParticle> *);

//! @endcond
