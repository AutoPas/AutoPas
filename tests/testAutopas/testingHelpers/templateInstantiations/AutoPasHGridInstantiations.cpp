/**
 * @file AutoPasHGridInstantiations.cpp
 *
 * HGrid-specific explicit template instantiations are kept in a separate translation unit
 * to reduce peak compile-time memory usage of AutoPasInstantiations.cpp.
 */

#include "autopas/AutoPasImpl.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/commonTypedefs.h"

//! @cond Doxygen_Suppress
template class autopas::AutoPas<HGridMockMolecule>;

template bool autopas::AutoPas<HGridMockMolecule>::computeInteractions(EmptyPairwiseFunctor<HGridMockMolecule> *);
//! @endcond
