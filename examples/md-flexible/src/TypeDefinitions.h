/**
 * @file TypeDefinitions.h
 * @author F. Gratl
 * @date 01.03.2021
 */

#pragma once

#if MD_FLEXIBLE_MODE == MULTISITE

#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
#include "molecularDynamicsLibrary/LJMultisiteFunctor.h"
#endif

#else

#include "molecularDynamicsLibrary/MoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
#include "molecularDynamicsLibrary/LJFunctor.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
#include "molecularDynamicsLibrary/LJFunctorAVX.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_HWY)
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_XSIMD)
#include "molecularDynamicsLibrary/LJFunctorXSIMD.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SIMDe)
#include "molecularDynamicsLibrary/LJFunctorSIMDe.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_MIPP)
#include "molecularDynamicsLibrary/LJFunctorMIPP.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SVE)
#include "molecularDynamicsLibrary/LJFunctorSVE.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC)
#include "molecularDynamicsLibrary/AxilrodTellerFunctor.h"
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

#if defined(MD_FLEXIBLE_FUNCTOR_HWY)
/**
 * Type of LJFunctorHWY used in md-flexible
 * Switches between mdLib::LJFunctorHWY and mdLib::LJMultisiteFunctorHWY as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 * @note mdLib::LJMultisiteFunctorHWY is yet to be written, so a compiler pre-processing error is thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have HWY support!"
#else

using LJFunctorTypeHWY = mdLib::LJFunctorHWY<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                             mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_XSIMD)
/**
 * Type of LJFunctorXSIMD used in md-flexible
 * Switches between mdLib::LJFunctorXSIMD and mdLib::LJMultisiteFunctorXSIMD as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 * @note mdLib::LJMultisiteFunctorXSIMD is yet to be written, so a compiler pre-processing error is thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have XSIMD support!"
#else

using LJFunctorTypeXSIMD = mdLib::LJFunctorXSIMD<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                                 mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_SIMDe)
/**
 * Type of LJFunctorSIMDe used in md-flexible
 * Switches between mdLib::LJFunctorSIMDe and mdLib::LJMultisiteFunctorSIMDe as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 * @note mdLib::LJMultisiteFunctorSIMDe is yet to be written, so a compiler pre-processing error is thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have SIMDe support!"
#else

using LJFunctorTypeSIMDe = mdLib::LJFunctorSIMDe<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                                 mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_MIPP)
/**
 * Type of LJFunctorMIPP used in md-flexible
 * Switches between mdLib::LJFunctorMIPP and mdLib::LJMultisiteFunctorMIPP as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 * @note mdLib::LJMultisiteFunctorMIPP is yet to be written, so a compiler pre-processing error is thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have MIPP support!"
#else

using LJFunctorTypeMIPP = mdLib::LJFunctorMIPP<ParticleType, true, true, autopas::FunctorN3Modes::Both,
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

#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC)
/**
 * Type of ATFunctor used in md-flexible.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "The Axilrod Teller functor does not have support for multisite molecules!"
#else
using ATFunctor = mdLib::AxilrodTellerFunctor<ParticleType, true, autopas::FunctorN3Modes::Both,
                                              mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
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
 * be used in the CMake.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#ifdef MD_FLEXIBLE_FUNCTOR_AUTOVEC
using LJFunctorTypeAbstract = mdLib::LJMultisiteFunctor<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS
using LJFunctorTypeAbstract = mdLib::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#endif

#else
#ifdef MD_FLEXIBLE_FUNCTOR_AUTOVEC
using LJFunctorTypeAbstract = mdLib::LJFunctor<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS
using LJFunctorTypeAbstract = mdLib::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#elif MD_FLEXIBLE_FUNCTOR_AVX
using LJFunctorTypeAbstract = mdLib::LJFunctorAVX<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_SVE
using LJFunctorTypeAbstract = mdLib::LJFunctorSVE<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_HWY
using LJFunctorTypeAbstract = mdLib::LJFunctorHWY<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_XSIMD
using LJFunctorTypeAbstract = mdLib::LJFunctorXSIMD<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_SIMDe
using LJFunctorTypeAbstract = mdLib::LJFunctorSIMDe<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_MIPP
using LJFunctorTypeAbstract = mdLib::LJFunctorMIPP<ParticleType, true, true>;
#endif

#endif
