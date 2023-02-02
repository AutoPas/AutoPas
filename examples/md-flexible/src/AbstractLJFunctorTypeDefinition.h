/**
* @file AbstractLJFunctorTypeDefinition.h
* @author S. Newcome
* @date 02/02/2023
*/

#pragma once

/**
 * We require access to a version of the force functor for non-simulation purposes, e.g. calculating FLOPs.
 * This is abstracted from whichever SoA implementation is used, so we pick any functor that is compiled.
 *
 * This file must be separate from TypeDefinitions.h, to avoid problems with compilations where LJFunctor Types cannot be
 * specified. In particular, for compiling mdFlexTests for the multi-site case, where there is only an AutoVec functor
 * implementation, which is by default disabled, and mdFlexTests has not got an option for enable this (and shouldn't).
*/
#ifdef MD_FLEXIBLE_USE_MULTI_SITE
#ifdef MD_FLEXIBLE_FUNCTOR_AUTOVEC
using LJFunctorTypeAbstract = autopas::LJMultisiteFunctor<ParticleType, true, true>;
#elif MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS
using LJFunctorTypeAbstract = autopas::LJMultisiteFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true>;
#else
#error "No Pairwise Force Functor has been compiled!"
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
#error "No Pairwise Force Functor has been compiled!"
#endif

#endif