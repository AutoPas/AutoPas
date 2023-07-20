/**
 * @file AutoPasClass.cpp
 *
 * Contains a explicit template instantiation for the main AutoPas class and the particle type used by md-flexible, as
 * determined by whether md-flexible is compiled with Multi-Site support or not. This is linked into the md-flexible
 * executable to enable the other compilation units to only declare, but not instantiate this template.
 */

#include "autopas/AutoPasImpl.h"
#include "src/TypeDefinitions.h"

template class autopas::AutoPas<ParticleType>;