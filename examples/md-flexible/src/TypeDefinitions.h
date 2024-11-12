/**
 * @file TypeDefinitions.h
 * @author F. Gratl
 * @date 01.03.2021
 */

#pragma once

#if MD_FLEXIBLE_MODE == MULTISITE

#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC) || defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
#include "molecularDynamicsLibrary/LJMultisiteFunctor.h"
#endif

#else

#include "molecularDynamicsLibrary/MoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC) || defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
#include "molecularDynamicsLibrary/LJFunctor.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
#include "molecularDynamicsLibrary/LJFunctorAVX.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SVE)
#include "molecularDynamicsLibrary/LJFunctorSVE.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_DEM)
#include "molecularDynamicsLibrary/GranularDEM.h"
#include "molecularDynamicsLibrary/DEMFunctor.h"
#endif

#endif

#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

/**
 * Precision used for particle representations. If you want to test other precisions change it here.
 */
using FloatPrecision = double;

/**
 * Type of the Particles used in md-flexible.
 * Switches between autopas::MoleculeLJ and autopas::MultisiteMoleculeLJ as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using ParticleType = mdLib::MultisiteMoleculeLJ;
#elif defined(MD_FLEXIBLE_FUNCTOR_DEM)
using ParticleType = mdLib::GranularDEM;
#else
using ParticleType = mdLib::MoleculeLJ;
#endif

namespace mdFlexibleTypeDefs {
/**
 * If AutoPas is compiled with FLOP logging enabled, use functors with FLOP counting enabled.
 */
constexpr bool countFLOPs =
#ifdef AUTOPAS_LOG_FLOPS
    true;
#else
    false;
#endif
}  // namespace mdFlexibleTypeDefs

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
/**
 * Type of LJFunctorTypeAutovec used in md-flexible.
 * Switches between mdLib::LJFunctor and mdLib::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAutovec = mdLib::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, false,
                                                       mdFlexibleTypeDefs::countFLOPs>;
#else
using LJFunctorTypeAutovec =
    mdLib::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, false, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
/**
 * Type of LJFunctorTypeAutovecGlobals used in md-flexible.
 * Switches between mdLib::LJFunctor and mdLib::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAutovecGlobals = mdLib::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                                              true, mdFlexibleTypeDefs::countFLOPs>;
#else
using LJFunctorTypeAutovecGlobals =
    mdLib::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
/**
 * Type of LJFunctorTypeAVX used in md-flexible.
 * Switches between mdLib::LJFunctorAVX and mdLib::LJMultisiteFunctorAVX as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 * @note mdLib::LJMultisiteFunctorAVX is yet to be written, so a compiler pre-processing error is thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have AVX support!"
#else
using LJFunctorTypeAVX =
    mdLib::LJFunctorAVX<ParticleType, true, true, autopas::FunctorN3Modes::Both, true, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SVE)
/**
 * Type of LJFunctorTypeSVE used in md-flexible.
 * Switches between mdLib::LJFunctorSVE and mdLib::LJMultisiteFunctorSVE as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 * @note mdLib::LJMultisiteFunctorSVE is yet to be written, so a compiler pre-processing error is thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have SVE support!"
#else
using LJFunctorTypeSVE =
    mdLib::LJFunctorSVE<ParticleType, false, true, autopas::FunctorN3Modes::Both, true, mdFlexibleTypeDefs::countFLOPs>; //TODO: set mixing, globalCalculations etc. to false for now
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_DEM)
/**
 * Type of LJFunctorTypeDEM used in md-flexible.
 */
#if MMD_FLEXIBLE_MODE == MULTISITE
#error "The DEM functor des not have support for multisite molecules!"
#else
using DEMFunctor = mdLib::DEMFunctor<ParticleType, true, autopas::FunctorN3Modes::Both, true, mdFlexibleTypeDefs::countFLOPs>;

#endif

#endif

/**
 * Type of the Particle Properties Library.
 * Set to the same precision as ParticleType.
 */
using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatPrecision, size_t>;
