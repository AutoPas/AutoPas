/**
 * @file loadEstimatorTest.cpp
 * @author fischerv
 * @date 04 May 2020
 */

#include "LoadEstimatorTest.h"

#include "autopas/containers/LoadEstimators.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

// using ::testing::_;

TEST(LoadEstimatorTest, testEqualDistributionSquaredParticlesPerCell) {
  std::array<size_t, 3> cellsPerDimension = {4, 4, 4};
  std::array<size_t, 3> particlesPerDimension = {8, 8, 8};
  std::array<double, 3> spacing = {.5, .5, .5};
  std::array<double, 3> offset = {.25, .25, .25};
  std::vector<FPCell> cells;
  cells.resize(cellsPerDimension[0] * cellsPerDimension[1] * cellsPerDimension[2]);

  autopasTools::generators::GridGenerator::fillWithParticles(cells, cellsPerDimension, particlesPerDimension,
                                                             typename FPCell::ParticleType(), spacing, offset);
  std::array<uint64_t, 3> lowerCorner = {0, 0, 0};

  for (uint64_t i = 1; i < 4; i++) {
    std::array<uint64_t, 3> upperCorner = {i, i, i};
    auto load = autopas::loadEstimators::squaredParticlesPerCell(cells, cellsPerDimension, lowerCorner, upperCorner);
    auto expectedLoad = (i + 1) * (i + 1) * (i + 1) * 64;
    EXPECT_EQ(load, expectedLoad);
  }
}

TEST(LoadEstimatorTest, testIncreasingDensitySquaredParticlesPerCell) {
  std::array<size_t, 3> cellsPerDimension = {3, 3, 3};
  std::vector<FPCell> cells;
  cells.resize(cellsPerDimension[0] * cellsPerDimension[1] * cellsPerDimension[2]);
  for (uint64_t i = 0; i < 3; i++) {
    std::array<size_t, 3> particlesPerDimension = {1, 3, 3 * (i + 1)};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0 / (i + 1)};
    std::array<double, 3> offset = {.5 + i, .5, 0.5 / (i + 1)};

    autopasTools::generators::GridGenerator::fillWithParticles(cells, cellsPerDimension, particlesPerDimension,
                                                               typename FPCell::ParticleType(), spacing, offset);
  }

  for (uint64_t i = 1; i < 3; i++) {
    std::array<uint64_t, 3> lowerCorner = {i, 0, 0};
    std::array<uint64_t, 3> upperCorner = {i, 2, 2};
    auto load = autopas::loadEstimators::squaredParticlesPerCell(cells, cellsPerDimension, lowerCorner, upperCorner);
    auto expectedLoad = 9 * (i + 1) * (i + 1);
    EXPECT_EQ(load, expectedLoad);
  }

  for (uint64_t i = 1; i < 3; i++) {
    std::array<uint64_t, 3> lowerCorner = {0, 0, i};
    std::array<uint64_t, 3> upperCorner = {2, 2, i};
    auto load = autopas::loadEstimators::squaredParticlesPerCell(cells, cellsPerDimension, lowerCorner, upperCorner);
    auto expectedLoad = 3 * (1 * 1 + 2 * 2 + 3 * 3);
    EXPECT_EQ(load, expectedLoad);
  }
}
