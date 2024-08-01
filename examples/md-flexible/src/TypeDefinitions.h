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

#if defined(MD_FLEXIBLE_FUNCTOR_AVX512_GS)
#include "molecularDynamicsLibrary/LJMultisiteFunctorAVX512_GS.h"
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX512_MASK)
#include "molecularDynamicsLibrary/LJMultisiteFunctorAVX512_Mask.h"
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

#include "molecularDynamicsLibrary/LJFunctorAVX512_Mask.h"

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

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
/**
 * Type of LJFunctorTypeAutovec used in md-flexible.
 * Switches between mdLib::LJFunctor and mdLib::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAutovec = mdLib::LJMultisiteFunctor<ParticleType, true>;
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
    mdLib::LJMultisiteFunctor<ParticleType, true, autopas::FunctorN3Modes::Both, true>;
#else
using LJFunctorTypeAutovecGlobals = mdLib::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX)
/**
 * Type of LJFunctorTypeAVX used in md-flexible.
 * Switches between mdLib::LJFunctorAVX and mdLib::LJMultisiteFunctorAVX512_GS as determined by CMake flag
 * MD_FLEXIBLE_MODE. The Multi-site variant uses cutoffs based on the distance between the center-of-masses.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have AVX support!. If your machine can use AVX512, consider using a variant of the AVX512 functor"
#else
using LJFunctorTypeAVX = mdLib::LJFunctorAVX<ParticleType, true, true>;
#endif

#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX512_GS)
/**
 * Type of LJFunctorTypeAVX512_GS used in md-flexible.
 * Switches between mdLib::LJFunctorAVX512_GS and mdLib::LJMultisiteFunctorAVX512_GS as determined by CMake flag MD_FLEXIBLE_MODE.
 * Handles cutoffs by using gather and scatter instruction to calculate forces only for interactions within the cutoff distance.
 * In the multi-site case, the cutoff checks are based on CoM-CoM distances.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAVX512_GS = mdLib::LJMultisiteFunctorAVX512_GS<ParticleType, false, autopas::FunctorN3Modes::Both, false, true>;
#else
#error "Single-Site Lennard-Jones Functor does not have AVX512 support!. "
#endif
#endif

#if defined(MD_FLEXIBLE_FUNCTOR_AVX512_MASK)
/**
 * Type of LJFunctorTypeAVX512_MASK used in md-flexible.
 * Switches between mdLib::LJFunctorAVX512_MASK and mdLib::LJMultisiteFunctorAVX512_MASK as determined by CMake flag MD_FLEXIBLE_MODE.
 * Handles cutoffs by using a mask for interactions beyond the cutoff distance. In the multi-site case, the cutoff checks
 * are based on site-to-site distances (i.e. every site pair is traversed).
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAVX512_MASK = mdLib::LJMultisiteFunctorAVX512_Mask<ParticleType, false, autopas::FunctorN3Modes::Both, false, true>;
#else
using LJFunctorTypeAVX512_MASK = mdLib::LJFunctorAVX512_Mask<ParticleType, true, true, autopas::FunctorN3Modes::Both, false, true>;
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
using LJFunctorTypeAbstract =
#if MD_FLEXIBLE_MODE == MULTISITE
#ifdef MD_FLEXIBLE_FUNCTOR_AUTOVEC
        mdLib::LJMultisiteFunctor<ParticleType, true>;
#elif MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS
        mdLib::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#elif MD_FLEXIBLE_FUNCTOR_AVX512_GS
        mdLib::LJMultisiteFunctorAVX512_GS<ParticleType, false, autopas::FunctorN3Modes::Both, false, true>;
#elif MD_FLEXIBLE_FUNCTOR_AVX512_MASK
        mdLib::LJMultisiteFunctorAV512_MASK<ParticleType, false, autopas::FunctorN3Modes::Both, false, true>;

#else
#include "molecularDynamicsLibrary/LJMultisiteFunctor.h"
    mdLib::LJMultisiteFunctor<ParticleType, true>;
#endif
#else
#ifdef MD_FLEXIBLE_FUNCTOR_AUTOVEC
        mdLib::LJFunctor<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS
        mdLib::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#elif MD_FLEXIBLE_FUNCTOR_AVX
        mdLib::LJFunctorAVX<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_SVE
        mdLib::LJFunctorSVE<ParticleType, true, true>;
#else
#include "molecularDynamicsLibrary/LJFunctor.h"
    mdLib::LJFunctor<ParticleType, true, true>;
#endif

#endif
