/**
 * @file commonTypedefs.h
 * @author F. Gratl
 * @date 6/21/18
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "autopas/particles/Particle.h"
#include "mocks/MockFunctor.h"

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
 * Short for the AutoPas single center Lennard-Jones molecule
 */
using Molecule = autopas::MoleculeLJ;
/**
 * Short for the Full Particle Cell with the single center Lennard-Jones molecule
 */
using FMCell = autopas::FullParticleCell<Molecule>;

// M prefix for mocks
/**
 * Short for Mock Functor
 */
using MFunctor = MockFunctor<autopas::Particle>;
