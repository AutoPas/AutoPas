/**
 * C08TraversalTest.h
 *
 *  Created on: 5/24/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <gtest/gtest.h>
#include <AutoPas.h>
#include <mocks/MockFunctor.h>

typedef MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> MFunctor;
typedef autopas::CellFunctor<autopas::Particle,
                             autopas::FullParticleCell<autopas::Particle>, MFunctor, false, true> MCellFunctor;
typedef autopas::FullParticleCell<autopas::Particle> FPCell;

class C08TraversalTest : public testing::Test {
 public:
  C08TraversalTest() = default;

  ~C08TraversalTest() = default;

  void fillWithParticles(std::vector<FPCell> &cells,
                         std::array<size_t, 3> particlesPerDim);
};
