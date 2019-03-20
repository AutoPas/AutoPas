/**
 * commonTypedefs.h
 *
 *  Created on: 6/21/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"
#include "mocks/MockFunctor.h"

// a place for typedefs that are commonly used in tests

typedef autopas::Particle Particle;
typedef autopas::FullParticleCell<autopas::Particle> FPCell;

typedef autopas::MoleculeLJ Molecule;
typedef autopas::FullParticleCell<autopas::MoleculeLJ> FMCell;

// M prefix for mocks
typedef MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> MFunctor;

typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor,
                             autopas::DataLayoutOption::aos, true>
    CellFunctorAoSN3;
typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor,
                             autopas::DataLayoutOption::aos, false>
    CellFunctorAoSNoN3;
typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor,
                             autopas::DataLayoutOption::soa, true>
    CellFunctorSoAN3;
typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor,
                             autopas::DataLayoutOption::soa, false>
    CellFunctorSoANoN3;
