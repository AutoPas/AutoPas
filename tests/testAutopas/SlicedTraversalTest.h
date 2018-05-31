/*
 * SlicedTraversalTest.h
 *
 *  Created on: 22 Jan 2018
 *      Author: gratl
 */

#pragma once

#include <AutoPas.h>
#include <gtest/gtest.h>
#include <mocks/MockFunctor.h>

typedef MockFunctor<autopas::Particle,
                    autopas::FullParticleCell<autopas::Particle>>
    MFunctor;
typedef autopas::CellFunctor<autopas::Particle,
                             autopas::FullParticleCell<autopas::Particle>,
                             MFunctor, false, true>
    MCellFunctor;
typedef autopas::FullParticleCell<autopas::Particle> FPCell;

class SlicedTraversalTest : public testing::Test {
 public:
  SlicedTraversalTest() = default;

  ~SlicedTraversalTest() override = default;

  void fillWithParticles(std::vector<FPCell> &cells,
                         std::array<size_t, 3> particlesPerDim);
};
