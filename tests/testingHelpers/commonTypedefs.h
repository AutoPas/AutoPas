/**
 * @file commonTypedefs.h
 * @author F. Gratl
 * @date 6/21/18
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "mocks/MockPairwiseFunctor.h"
#include "mocks/MockTriwiseFunctor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

// a place for aliases that are commonly used in tests

/**
 * Short for AutoPas ParticleBaseFP64
 */
using ParticleFP64 = autopas::ParticleBaseFP64;
/**
 * Short for a FullParticle Cell with the AutoPas ParticleFP64
 */
using FPCell = autopas::FullParticleCell<ParticleFP64>;

/**
 * Short for the AutoPas single site Lennard-Jones molecule
 */
using Molecule = mdLib::MoleculeLJ;
/**
 * Short for the Full Particle Cell with the single center Lennard-Jones molecule
 */
using FMCell = autopas::FullParticleCell<Molecule>;

// M prefix for mocks
/**
 * Short for Mock Pairwise Functor
 */
using MPairwiseFunctor = MockPairwiseFunctor<ParticleFP64>;
/**
 * Short for Mock Triwise Functor
 */
using MTriwiseFunctor = MockTriwiseFunctor<ParticleFP64>;

/**
 * Helper alias for LJFunctor, with more defaults geared towards testing.
 * This facilitates writing tests and tries to reduce the number of template instantiations.
 */
template <bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
using LJFunctorType =
    mdLib::LJFunctor<Molecule, applyShift, useMixing, useNewton3, calculateGlobals, countFLOPs, relevantForTuning>;

/**
 * Helper alias for specialization of LJFunctorType with globals and shift enabled but mixing disabled.
 */
using LJFunctorGlobals = LJFunctorType</* shifting */ true, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                                       /*globals*/ true>;
