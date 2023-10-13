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
#elif defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
#include "molecularDynamicsLibrary/LJAbsoluteMultiSiteFunctor.h"
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
#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
using ParticleType = mdLib::AbsoluteMultiSiteMoleculeLJ;
#else
using ParticleType = mdLib::MultisiteMoleculeLJ;
#endif
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
using LJFunctorTypeAutovec = mdLib::LJFunctor<ParticleType, true, true>;
#endif

#endif

//@todo (johnny) considering we are defining LJFunctorTypeAbstract further down i don't even know why we are doing this, i am just pattern matching
#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
/**
* Type of LJFunctorTypeAutovec used in md-flexible.
* Can only be LJAbsoluteMultiSiteFunctor for MultiSite molecules. A corresponding functor for single site molecules wouldn't really make sense
 */
#if MD_FLEXIBLE_MODE == MULTISITE
using LJFunctorTypeAbsPos = mdLib::LJAbsoluteMultiSiteFunctor<ParticleType, true, true>;
#else
#error "Single site simulation can't use absolute position multisite functor."
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
#elif MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS
using LJFunctorTypeAbstract = mdLib::LJAbsoluteMultiSiteFunctor<ParticleType, true, true>;
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
