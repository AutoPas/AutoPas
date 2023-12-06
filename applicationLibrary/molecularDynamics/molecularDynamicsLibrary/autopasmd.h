/**
 * @file autopasmd.h
 * Main include file for all molecular dynamics functionality for convenience.
 *
 * @author seckler
 * @date 15.01.20
 */

#pragma once

#include "LJFunctor.h"
#include "MieFunctor.h"
#ifdef __AVX__
#include "LJFunctorAVX.h"
#include "MieFunctorAVX.h"
#include "MieFunctorAVXFixedExponents.h"
#endif
#ifdef __ARM_FEATURE_SVE
#include "LJFunctorSVE.h"
#include "MieFunctorSVE.h"
#endif
#include "MoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
