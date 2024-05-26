/**
 * @file AutoPasInstantiations.cpp
 *
 * Contains a explicit template instantiation so that they are not generated on a test-by-test basis.
 * Generating similar templates in the same translation unit also saves memory.
 */

#include "autopas/AutoPasImpl.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/EmptyFunctor.h"
#include "testingHelpers/NonConstructibleParticle.h"
#include "testingHelpers/commonTypedefs.h"

//! @cond Doxygen_Suppress
template class autopas::AutoPas<Molecule>;
template class autopas::AutoPas<NonConstructibleParticle>;

template bool autopas::AutoPas<Molecule>::iteratePairwise(mdLib::LJFunctor<Molecule> *);
template bool autopas::AutoPas<Molecule>::iteratePairwise(
    LJFunctorType</* shifting */ true, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                  /*globals*/ true> *);
template bool autopas::AutoPas<Molecule>::iteratePairwise(
    LJFunctorType</* shifting */ true, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                  /*globals*/ false> *);
template bool autopas::AutoPas<Molecule>::iteratePairwise(EmptyFunctor<Molecule> *);
template bool autopas::AutoPas<NonConstructibleParticle>::iteratePairwise(MockFunctor<NonConstructibleParticle> *);

//! @endcond
