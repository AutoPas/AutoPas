/**
 * @file commonTypedefs.h
 * @author F. Gratl
 * @date 6/21/18
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "mocks/MockPairwiseFunctor.h"
#include "mocks/MockTriwiseFunctor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

// a place for usings that are commonly used in tests

/**
 * Short for AutoPas Particle
 */
using Particle = autopas::Particle;
/**
 * Short for a FullParticle Cell with the AutoPas Particle
 */
using FPCell = autopas::FullParticleCell<autopas::Particle>;

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
using MPairwiseFunctor = MockPairwiseFunctor<autopas::Particle>;
/**
 * Short for Mock Triwise Functor
 */
using MTriwiseFunctor = MockTriwiseFunctor<autopas::Particle>;

namespace autopasTestingTypeDefs {
/**
 * If AutoPas is compiled with FLOP logging enabled, use functors with FLOP counting enabled.
 */
constexpr bool countFLOPs =
#ifdef AUTOPAS_LOG_FLOPS
    true;
#else
    false;
#endif
}  // namespace autopasTestingTypeDefs

/**
 * Helper alias for LJFunctor, which specifies Particle as Molecule (as defined above) and countFLOPs as defined above.
 * Primarily, this is useful to avoid specifying out default template when we just want to change countFLOPs, which we
 * need to do to avoid warnings/exceptions when the tests are compiled with AUTOPAS_LOG_FLOPS=ON.
 */
template <bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
using LJFunctorType = mdLib::LJFunctor<Molecule, applyShift, useMixing, useNewton3, calculateGlobals,
                                       autopasTestingTypeDefs::countFLOPs, relevantForTuning>;

/**
 * Helper alias for specialization of LJFunctorType with globals and shift enabled but mixing disabled.
 */
using LJFunctorGlobals = LJFunctorType</* shifting */ true, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                                       /*globals*/ true>;
