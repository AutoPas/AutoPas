/**
 * @file AutoPasMolecule.cpp
 *
 * Contains a explicit template instantiation for the main AutoPas class and the Molecule data type. This
 * is linked into the test executable to enable the other compilation units to only declare, but not instantiate
 * this template.
 */

#include "autopas/AutoPasImpl.h"
#include "testingHelpers/NonConstructibleParticle.h"

template class autopas::AutoPas<NonConstructibleParticle>;