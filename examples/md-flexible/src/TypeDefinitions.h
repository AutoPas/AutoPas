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

#if defined(MD_FLEXIBLE_FUNCTOR_AT)
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

#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
/**
 * Type of LJFunctorTypeAutovec used in md-flexible.
 * Switches between mdLib::LJFunctor and mdLib::LJMultisiteFunctor as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAutovec = mdLib::LJMultisiteFunctor<ParticleType, true, true>;
#else
using LJFunctorTypeAutovec = mdLib::LJFunctor<ParticleType, true, true, false>; // Changed from true
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
 * Switches between mdLib::LJFunctorAVX and mdLib::LJMultisiteFunctorAVX as determined by CMake flag
 * MD_FLEXIBLE_MODE.
 * @note mdLib::LJMultisiteFunctorAVX is yet to be written, so a compiler pre-processing error is thrown.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "Multi-Site Lennard-Jones Functor does not have AVX support!"
#else
using LJFunctorTypeAVX = mdLib::LJFunctorAVX<ParticleType, true, true>;
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

#if defined(MD_FLEXIBLE_FUNCTOR_AT)
/**
 * Type of LJFunctorTypeAT used in md-flexible.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
#error "The Axilrod Teller functor does not have support for multisite molecules!"
#else
using ATFunctor = mdLib::AxilrodTellerFunctor<ParticleType, true, true>; //Changed
#endif

#endif


/**
 * Look-up Table types
 */

#include "molecularDynamicsLibrary/LJLookUpTable.h"

#if defined(MD_USE_LJ_LUT)

const bool USE_LJ_LUT = true;

#if LJ_LUT_INTERVALL == EVEN_SPACING

#if LJ_LUT_INTERPOLATION_FUNCTION == NEXT_NEIGHBOR

using LJLookUpTableType = ForceLookUpTable::LJLookUpTable<ForceLookUpTable::evenSpacing, ForceLookUpTable::nextNeighbor, FloatPrecision, size_t>;

#elif LJ_LUT_INTERPOLATION_FUNCTION

using LJLookUpTableType = ForceLookUpTable::LJLookUpTable<ForceLookUpTable::evenSpacing, ForceLookUpTable::linear, FloatPrecision, size_t>;

#endif // LJ_LUT_INTERPOLATION_FUNCTION

#endif // LJ_LUT_INTERVALL

#else

using LJLookUpTableType = ForceLookUpTable::LJLookUpTable<ForceLookUpTable::evenSpacing, ForceLookUpTable::nextNeighbor, FloatPrecision, size_t>;

const bool USE_LJ_LUT =  false;

#endif // MD_USE_LJ_LUT

#include "molecularDynamicsLibrary/ATLookUpTable.h"

#if defined(MD_USE_AT_LUT)

const bool USE_AT_LUT = true;

#if AT_LUT_INTERVALL == EVEN_SPACING

#if AT_LUT_INTERPOLATION_FUNCTION == NEXT_NEIGHBOR

using ATLookUpTableType = ForceLookUpTable::ATLookUpTable<ForceLookUpTable::evenSpacing, ForceLookUpTable::nextNeighbor, FloatPrecision, size_t>;

#elif AT_LUT_INTERPOLATION_FUNCTION == LINEAR

using ATLookUpTableType = ForceLookUpTable::ATLookUpTable<ForceLookUpTable::evenSpacing, ForceLookUpTable::linear, FloatPrecision, size_t;

#endif // AT_LUT_INTERPOLATION_FUNCTION

#endif // AT_LUT_INTERVALL
Oh,
#else

using ATLookUpTableType = ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing, ForceLookUpTable::nextNeighbor, FloatPrecision, size_t>;

const bool USE_AT_LUT = false;

#endif // MD_USE_AT_LUT

/**
 * Type of the Particle Properties Library.
 * Set to the same precision as ParticleType.
 */
using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatPrecision, size_t>;

/**
 * We require access to a version of the force functor for non-computeInteractions purposes, e.g. calculating FLOPs or
 * AoS functor calls. This is abstracted from whichever SoA implementation is used, so we pick any functor that is
 * chosen to be used in the CMake.
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

#ifdef MD_FLEXIBLE_FUNCTOR_AT
using ATFunctorTypeAbstract = mdLib::AxilrodTellerFunctor<ParticleType, true, true>; // Changed
#endif

#endif
