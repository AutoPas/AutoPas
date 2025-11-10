/**
 * @file TypeDefinitions.h
 * @author F. Gratl
 * @date 01.03.2021
 */

#pragma once
#include "molecularDynamicsLibrary/AxilrodTellerMutoMultisiteFunctor.h"
#include "molecularDynamicsLibrary/MethaneMultisitePairwiseFunctor.h"

#if MD_FLEXIBLE_MODE == MULTISITE

#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
#include "molecularDynamicsLibrary/LJMultisiteFunctor.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_METHANE)
#include "molecularDynamicsLibrary/MethaneMultisitePairwiseFunctor.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC)
#include "molecularDynamicsLibrary/AxilrodTellerMutoMultisiteFunctor.h"
#endif

#else

#include "molecularDynamicsLibrary/MoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
#include "molecularDynamicsLibrary/LJFunctor.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
#include "molecularDynamicsLibrary/LJFunctorAVX.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SVE)
#include "molecularDynamicsLibrary/LJFunctorSVE.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC)
#include "molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h"
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

/**
 * If md-flexible is compiled with globals calculations enabled, use functors which calculate globals.
 */
constexpr bool calcGlobals =
#ifdef MD_FLEXIBLE_CALC_GLOBALS
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
using LJFunctorTypeAutovec = mdLib::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                                       mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
#else
using LJFunctorTypeAutovec = mdLib::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                              mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
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
using LJFunctorTypeAVX = mdLib::LJFunctorAVX<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                             mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
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
using LJFunctorTypeSVE = mdLib::LJFunctorSVE<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                             mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC)
/**
 * Type of ATMFunctor used in md-flexible.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using ATMFunctor =
    mdLib::AxilrodTellerMutoMultisiteFunctor<ParticleType, true, autopas::FunctorN3Modes::Both,
                                             mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
#else
using ATMFunctor = mdLib::AxilrodTellerMutoFunctor<ParticleType, true, autopas::FunctorN3Modes::Both,
                                                   mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

// Methane pairwise functor
#if defined(MD_FLEXIBLE_FUNCTOR_METHANE)
#if MD_FLEXIBLE_MODE == MULTISITE
using MethaneFunctor =
    mdLib::MethaneMultisitePairwiseFunctor<ParticleType, autopas::FunctorN3Modes::Both, mdFlexibleTypeDefs::calcGlobals,
                                           mdFlexibleTypeDefs::countFLOPs>;
#else
#error \
    "There exists no pairwise functor for single-site methane molecules! Consider using LJ or multisite mode instead."
#endif
#endif

/**
 * Type of the Particle Properties Library.
 * Set to the same precision as ParticleType.
 */
using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatPrecision, size_t>;
