/**
 * @file SimulationClass.cpp
 *
 * Contains a explicit template instantiation for the main Simulation class and the particle types used by md-flexible. This
 * is linked into the md-flexible executable to enable the other compilation units to only declare, but not instantiate
 * this template.
 */

#include "src/Simulation.h"
#include "src/TypeDefinitions.h"

template class Simulation<SingleSiteMolecule>;
template class Simulation<MultiSiteMolecule>;