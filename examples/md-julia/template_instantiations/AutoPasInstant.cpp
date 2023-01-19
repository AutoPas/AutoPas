/**
 * @file AutoPasInstant.cpp
 *
 * Contains a explicit template instantiation for the main AutoPas class and the particle type used by md-flexible. This
 * is linked into the md-flexible executable to enable the other compilation units to only declare, but not instantiate
 * this template.
 */
#pragma once
#include "autopas/AutoPasImpl.h"
#include "src/TypeDef.h"

// template class autopas::AutoPas<ParticleType>;