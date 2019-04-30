/**
 * @file commonTypedefs.h
 * @author F. Gratl
 * @date 6/21/18
 */

#pragma once

#include "mocks/MockFunctor.h"

// a place for typedefs that are commonly used in tests

/**
 * Short for AutoPas Particle
 */
typedef autopas::Particle Particle;
/**
 * Short for a FullParticle Cell with the AutoPas Particle
 */
typedef autopas::FullParticleCell<autopas::Particle> FPCell;

/**
 * Short for the AutoPas single center Lennard-Jones molecule
 */
typedef autopas::MoleculeLJ Molecule;
/**
 * Short for the Full Particle Cell with the single center Lennard-Jones molecule
 */
typedef autopas::FullParticleCell<autopas::MoleculeLJ> FMCell;

// M prefix for mocks
/**
 * Short for Mock Functor
 */
typedef MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> MFunctor;

/**
 * Short for the Cell Functor with AoS and Newton 3
 */
typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor, false, true>
    CellFunctorAoSN3;
/**
 * Short for the Cell Functor with AoS and no Newton 3
 */
typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor, false, false>
    CellFunctorAoSNoN3;
/**
 * Short for the Cell Functor with SoA and Newton 3
 */
typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor, true, true>
    CellFunctorSoAN3;
/**
 * Short for the Cell Functor with SoA and no Newton 3
 */
typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor, true, false>
    CellFunctorSoANoN3;
