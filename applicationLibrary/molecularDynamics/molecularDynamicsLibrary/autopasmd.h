/**
 * @file autopasmd.h
 * Main include file for all molecular dynamics functionality for convenience.
 *
 * @author seckler
 * @date 15.01.20
 */

#pragma once

#include "LJFunctor.h"
#ifdef __AVX__
#include "LJFunctorAVX.h"
#endif
#ifdef __ARM_FEATURE_SVE
#include "LJFunctorSVE.h"
#endif
#include "MoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
