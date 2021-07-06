/**
 * @file OctreeTests.h
 * @author Johannes Spies
 * @date 15.04.2021
 */

#pragma once

#include <gtest/gtest.h>
#include <testingHelpers/commonTypedefs.h>

#include <array>
#include <vector>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "mocks/MockFunctor.h"

/**
 * A pair of particle counts and halo particle counts.
 */
using GeneratorSpec =
    std::tuple<std::array<double, 3> /*boxMax*/, int unsigned /*numParticles*/, int unsigned /*numHaloParticles*/>;

class OctreeTest : public AutoPasTestBase, public ::testing::WithParamInterface<GeneratorSpec> {
 public:
  using Vector3DList = std::vector<std::array<double, 3>>;

  std::pair<Vector3DList, std::vector<std::tuple<unsigned long, unsigned long, double>>> calculateForcesAndPairs(
      autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
      autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, size_t numParticles,
      size_t numHaloParticles, std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cellSizeFactor,
      double cutoff, double skin, double interactionLength, Vector3DList particlePositions,
      Vector3DList haloParticlePositions);

 public:
  MockFunctor<Molecule> mockFunctor;
};
