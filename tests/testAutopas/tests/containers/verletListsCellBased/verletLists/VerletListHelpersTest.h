/**
 * @file VerletListHelpersTest.cpp
 * @author muehlhaeusser
 * @date 2026-07-26
 */

#pragma once

#include <gtest/gtest.h>

#include <unordered_map>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/particles/ParticleDefinitions.h"

using namespace autopas;

class VerletListHelpersTest : public ::testing::Test {
 protected:
  using ParticleType = ParticleBaseFP64;
  using CellType = FullParticleCell<ParticleType>;
  using Helpers = VerletListHelpers<ParticleType>;

  void SetUp() override {}

  /**
   * Helper function to generate N interacting particles at the origin.
   */
  void generateDenseParticles(size_t n) {
    for (size_t i = 0; i < n; ++i) {
      // All particles at origin, so distance is 0.0 (they all interact)
      cell.addParticle(ParticleType({0., 0., 0.}, {0., 0., 0.}, i));
      particleToIndex[&cell[i]] = i;
    }
  }

  CellType cell;
  std::unordered_map<const ParticleType *, size_t> particleToIndex;
  const double interactionLength = 1.0;
};
