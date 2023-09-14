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

#if defined(MD_FLEXIBLE_FUNCTOR_AVX) || defined(MD_FLEXIBLE_FUNCTOR_AVX_GS)
#include "molecularDynamicsLibrary/LJMultisiteFunctorAVX.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX_STS)
#include "molecularDynamicsLibrary/LJMultisiteFunctorAVX_STS.h"
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

/**
 * Some Functors require a vecLength i.e. the number of doubles which fit into a floating register.
 * Currently, this is by default set to 4. If AVX512 is used this is set to 8.
 * Todo check SVE
 */
const size_t vecLength =
#ifdef __AVX512__
    8;
#else
    4;
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
/**
 * Type of LJFunctorTypeAutovec used in md-flexible.
 * Switches between mdLib::LJFunctor and mdLib::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAutovec = mdLib::LJMultisiteFunctor<ParticleType, true, true>;
#else
using LJFunctorTypeAutovec = mdLib::LJFunctor<ParticleType, true, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
/**
 * Type of LJFunctorTypeAutovecGlobals used in md-flexible.
 * Switches between mdLib::LJFunctor and mdLib::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAutovecGlobals =
    mdLib::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#else
using LJFunctorTypeAutovecGlobals = mdLib::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
/**
 * Type of LJFunctorTypeAVX used in md-flexible.
 * Switches between mdLib::LJFunctorAVX and mdLib::LJMultisiteFunctorAVX (with 0/1-Masks) as determined by CMake flag
 * MD_FLEXIBLE_MODE. The Multi-site variant uses cutoffs based on the distance between the center-of-masses.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAVX = mdLib::LJMultisiteFunctorAVX<ParticleType, false, true, autopas::FunctorN3Modes::Both, false, true, true, vecLength>;
#else
using LJFunctorTypeAVX = mdLib::LJFunctorAVX<ParticleType, true, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX_GS)
/**
 * Type of LJFunctorTypeAVXGS used in md-flexible.
 * Switches between a single-site LJ Functor with Gather/Scatter instructions and mdLib::LJMultisiteFunctorAVX (with G/S) as determined by CMake flag MD_FLEXIBLE_MODE.
 * The Multi-site variant uses cutoffs based on the distance between the center-of-masses.
 * @note No Single-site LJ Functor using Gather/Scatter instructions has been created, so a preprocessor error will be thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAVXGS = mdLib::LJMultisiteFunctorAVX<ParticleType, false, true, autopas::FunctorN3Modes::Both, false, true, false, vecLength>;
#else
#error "Gather/Scatter LJ Functor does not have single-site support!"
#endif
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX_STS)
/**
 * Type of LJFunctorTypeAVXSTS used in md-flexible.
 * The only difference between this and LJFunctorTypeAVX is the cutoff criterion - which here is based individually on
 * site to site distances. There is no distinction between this functor and LJFunctorTypeAVX in single-site mode, and
 * thus a preprocessor error is thrown if used in that mode.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAVXGSSTS = mdlib::LJMultisiteFunctorAVX_STS<ParticleType, false, true, autopas::FunctorN3Modes::Both, false, true, true, vecLength>;
#else
#error "MD_FLEXIBLE_FUNCTOR_AVX_STS is not valid in SINGLESITE mode!"
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
using LJFunctorTypeSVE = mdLib::LJFunctorSVE<ParticleType, true, true>;
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
#endif

#endif
