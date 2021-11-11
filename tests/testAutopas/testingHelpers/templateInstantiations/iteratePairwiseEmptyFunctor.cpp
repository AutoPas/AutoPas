/**
 * @file iteratePairwiseLJFunctor.cpp
 *
 * Contains a explicit template instantiation for the iteratePairwise() method with the EmptyFunctor of the main
 * AutoPas class and a particle type. This is linked into the test executable to enable the
 * other compilation units to only declare, but not instantiate this template.
 */

#include "autopas/AutoPasImpl.h"
#include "testingHelpers/EmptyFunctor.h"
#include "testingHelpers/commonTypedefs.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<Molecule>::iteratePairwise(EmptyFunctor<Molecule> *);
//! @endcond
