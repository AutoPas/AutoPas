/**
 * @file TypeDefinitions.h
 * @author F. Gratl
 * @date 01.03.2021
 */

#pragma once

#if defined(MD_FLEXIBLE_USE_MULTI_SITE)

#include "autopas/molecularDynamics/MultisiteMoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC) || defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
#include "autopas/molecularDynamics/LJMultisiteFunctor.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
//#include "autopas/molecularDynamics/LJMultisiteFunctorCTS.h"
#include "autopas/molecularDynamics/LJMultisiteFunctorAVX.h"
#endif

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
 * Switches between autopas::MoleculeLJ and autopas::MultisiteMoleculeLJ as determined by CMake flag
 * MD_FLEXIBLE_USE_MULTI_SITE.
 */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
using ParticleType = autopas::MultisiteMoleculeLJ;
#else
using ParticleType = autopas::MoleculeLJ;
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
/**
 * Type of LJFunctorTypeAutovec used in md-flexible.
 * Switches between autopas::LJFunctor and autopas::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_USE_MULTI_SITE.
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
 * Switches between autopas::LJFunctor and autopas::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_USE_MULTI_SITE.
 */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
using LJFunctorTypeAutovecGlobals =
    autopas::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#else
using LJFunctorTypeAutovecGlobals = autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
/**
 * Type of LJFunctorTypeAVX used in md-flexible.
 * Switches between autopas::LJFunctorAVX and autopas::LJMultisiteFunctorAVX as determined by CMake flag
 * MD_FLEXIBLE_USE_MULTI_SITE.
 * @note autopas::LJMultisiteFunctorAVX is yet to be written, so a compiler pre-processing error is thrown.
 */
#if defined(MD_FLEXIBLE_USE_MULTI_SITE)
//#pragma message "Multi-Site Lennard-Jones Functor with AVX is currently WIP!"
using LJFunctorTypeAVX = autopas::LJMultisiteFunctorAVX<ParticleType, true, true>;
//using LJFunctorTypeAVX = autopas::LJMultisiteFunctorCTS<ParticleType, true, true>;
#else
using LJFunctorTypeAVX = autopas::LJFunctorAVX<ParticleType, true, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SVE)
/**
 * Type of LJFunctorTypeSVE used in md-flexible.
 * Switches between autopas::LJFunctorSVE and autopas::LJMultisiteFunctorSVE as determined by CMake flag
 * MD_FLEXIBLE_USE_MULTI_SITE.
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

/**
 * We require access to a version of the force functor for non-iteratePairwise purposes, e.g. calculating FLOPs or AoS
 * functor calls. This is abstracted from whichever SoA implementation is used, so we pick any functor that is chosen to
 * be used in the CMake. If no (valid) implementation is chosen, this is set to some arbitrary valid implementation,
 * e.g. AutoVec.
 */
#ifdef MD_FLEXIBLE_USE_MULTI_SITE
#ifdef MD_FLEXIBLE_FUNCTOR_AUTOVEC
using LJFunctorTypeAbstract = autopas::LJMultisiteFunctor<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS
using LJFunctorTypeAbstract =
    autopas::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#else
#include "autopas/molecularDynamics/LJMultisiteFunctor.h"
using LJFunctorTypeAbstract = autopas::LJMultisiteFunctor<ParticleType, true, true>;
#endif

#else
#ifdef MD_FLEXIBLE_FUNCTOR_AUTOVEC
using LJFunctorTypeAbstract = autopas::LJFunctor<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS
using LJFunctorTypeAbstract = autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#elif MD_FLEXIBLE_FUNCTOR_AVX
using LJFunctorTypeAbstract = autopas::LJFunctorAVX<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_SVE
using LJFunctorTypeAbstract = autopas::LJFunctorSVE<ParticleType, true, true>;
#else
#include "autopas/molecularDynamics/LJFunctor.h"
using LJFunctorTypeAbstract = autopas::LJFunctor<ParticleType, true, true>;
#endif

#endif
