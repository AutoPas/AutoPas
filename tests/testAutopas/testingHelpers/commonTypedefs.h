/**
 * commonTypedefs.h
 *
 *  Created on: 6/21/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include "mocks/MockFunctor.h"

// a place for typedefs that are commonly used in tests

typedef autopas::Particle Particle;

typedef autopas::FullParticleCell<autopas::Particle> FPCell;

typedef MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> MFunctor;

typedef autopas::CellFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>, MFunctor, false, true>
    MCellFunctorAoSN3;

