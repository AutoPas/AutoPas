/**
 * @file AutoPasClass.cpp
 *
 * Contains a explicit template instantiation for the main AutoPas class and the particle types used by md-flexible. This
 * is linked into the md-flexible executable to enable the other compilation units to only declare, but not instantiate
 * this template.
 */

#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

template class autopas::AutoPas<SingleSiteMolecule>;
template class autopas::AutoPas<MultiSiteMolecule>;