/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 01.03.2021
 */

#pragma once

#if defined(MD_FLEXIBLE_USE_MULTI_SITE)

#include "autopas/molecularDynamics/MultisiteMoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC) || defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
#include "autopas/molecularDynamics/LJMultisiteFunctor.h"
#endif

#else

#include "autopas/molecularDynamics/MoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC) || defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
#include "autopas/molecularDynamics/LJFunctor.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SVE)
#include "autopas/molecularDynamics/LJFunctorSVE.h"
#endif

#endif

#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"

/**
 * Precision used for particle representations. If you want to test other precisions change it here.
 */
using FloatPrecision = double;

/**
 * Type of the Particles used in md-flexible.
 * Switches between autopas::MoleculeLJ and autopas::MultisiteMoleculeLJ as determined by CMake flag MD_FLEXIBLE_USE_MULTI_SITE.
 */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
  using ParticleType = autopas::MultisiteMoleculeLJ;
#else
  using ParticleType = autopas::MoleculeLJ;
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
/**
 * Type of LJFunctorTypeAutovec used in md-flexible.
 * Switches between autopas::LJFunctor and autopas::LJMultisiteFunctor as determined by CMake flag MD_FLEXIBLE_USE_MULTI_SITE.
 */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
  using LJFunctorTypeAutovec = autopas::LJMultisiteFunctor<ParticleType, true, true>;
#else
  using LJFunctorTypeAutovec = autopas::LJFunctor<ParticleType, true, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
/**
 * Type of LJFunctorTypeAutovecGlobals used in md-flexible.
 * Switches between autopas::LJFunctor and autopas::LJMultisiteFunctor as determined by CMake flag MD_FLEXIBLE_USE_MULTI_SITE.
 */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
  using LJFunctorTypeAutovecGlobals = autopas::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#else
  using LJFunctorTypeAutovecGlobals = autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
/**
 * Type of LJFunctorTypeAVX used in md-flexible.
 * Switches between autopas::LJFunctorAVX and autopas::LJMultisiteFunctorAVX as determined by CMake flag MD_FLEXIBLE_USE_MULTI_SITE.
 * @note autopas::LJMultisiteFunctorAVX is yet to be written, so a compiler pre-processing error is thrown.
 */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
#error "Multi-Site Lennard-Jones Functor does not have AVX support!"
#else
  using LJFunctorTypeAVX = autopas::LJFunctorAVX<ParticleType, true, true>;
#endif

#endif


#if defined(MD_FLEXIBLE_FUNCTOR_SVE)
  /**
 * Type of LJFunctorTypeSVE used in md-flexible.
 * Switches between autopas::LJFunctorSVE and autopas::LJMultisiteFunctorSVE as determined by CMake flag MD_FLEXIBLE_USE_MULTI_SITE.
 * @note autopas::LJMultisiteFunctorSVE is yet to be written, so a compiler pre-processing error is thrown.
   */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
#error "Multi-Site Lennard-Jones Functor does not have SVE support!"
#else
  using LJFunctorTypeSVE = autopas::LJFunctorSVE<ParticleType, true, true>;
#endif

#endif

/**
 * Type of the Particle Properties Library.
 * Set to the same precision as ParticleType.
 */
using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatPrecision, size_t>;
