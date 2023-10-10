/**
 * @file autopasmd.h
 * Main include file for all molecular dynamics functionality for convenience.
 *
 * @author seckler
 * @date 15.01.20
 */

#pragma once

#include "autopas/molecularDynamics/LJFunctor.h"
#ifdef __AVX__
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#endif
#ifdef __ARM_FEATURE_SVE
#include "autopas/molecularDynamics/LJFunctorSVE.h"
#endif
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
