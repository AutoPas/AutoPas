/**
 * @file MoleculeBase.h
 * @date 02 Oct 2026
 * @author F. Duhr
 */

#pragma once

#include "autopas/particles/ParticleDefinitions.h"

namespace mdLib {

/**
 * Base class of all molecules. The floating point precision is chosen via the CMake option MD_USE_FLOAT_PRECISION:
 * float if it is set, double otherwise.
 */
#ifdef MD_USE_FLOAT_PRECISION
using MoleculeBase = autopas::ParticleBaseFP32;
#else
using MoleculeBase = autopas::ParticleBaseFP64;
#endif

}  // namespace mdLib
