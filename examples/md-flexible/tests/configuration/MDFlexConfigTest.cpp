/**
 * @file MDFlexConfigTest.cpp
 * @author F. Gratl
 * @date 04/06/2020
 */

#include "MDFlexConfigTest.h"

#include <gmock/gmock-matchers.h>

using ::testing::ElementsAreArray;

TEST_F(MDFlexConfigTest, GridBoxMinMax) {
  double spacing = 1.2;
  std::array<size_t, 3> particlesPerDim{4, 5, 6};
  std::array<double, 3> bottomLeftCorner{0, 1, 2};
  CubeGrid cubeGrid(zero, 0, 1, 1, 1, particlesPerDim, spacing, bottomLeftCorner);

  _configuration->cubeGridObjects.push_back(cubeGrid);
  _configuration->calcSimulationBox();

  EXPECT_THAT(_configuration->boxMin.value, ElementsAreArray({-spacing / 2, 0., 0.}));

  std::array<double, 3> expectedBoxMax;
  for (unsigned int i = 0; i < expectedBoxMax.size(); ++i) {
    // position of last particle
    expectedBoxMax[i] = (particlesPerDim[i] - 1) * spacing + bottomLeftCorner[i];
    // + offset so that object is not exactly at the edge of the domain (important for periodics)
    expectedBoxMax[i] += spacing * 0.5;
  }
  EXPECT_THAT(_configuration->boxMax.value, ElementsAreArray(expectedBoxMax));
}

TEST_F(MDFlexConfigTest, SphereBoxMinMax) {
  double spacing = 1.2;
  size_t radius = 5;
  std::array<double, 3> center{0, 1, 2};
  Sphere sphere(zero, 0, 1, 1, 1, center, radius, spacing);

  _configuration->sphereObjects.push_back(sphere);
  _configuration->calcSimulationBox();

  std::array<double, 3> expectedBoxMax;
  std::array<double, 3> expectedBoxMin;
  for (unsigned int i = 0; i < expectedBoxMax.size(); ++i) {
    // position of last particle
    expectedBoxMax[i] = center[i] + radius * spacing;
    expectedBoxMin[i] = center[i] - radius * spacing;
    // + offset so that object is not exactly at the edge of the domain (important for periodics)
    expectedBoxMax[i] += spacing * 0.5;
    expectedBoxMin[i] -= spacing * 0.5;
  }
  EXPECT_THAT(_configuration->boxMin.value, ElementsAreArray(expectedBoxMin));
  EXPECT_THAT(_configuration->boxMax.value, ElementsAreArray(expectedBoxMax));
}

TEST_F(MDFlexConfigTest, GaussBoxMinMax) {
  size_t numParticles = 10;
  std::array<double, 3> boxLength{4, 5, 6};
  std::array<double, 3> bottomLeftCorner{-1, 0, 1};
  CubeGauss cubeGauss(zero, 0, 1, 1, 1, numParticles, boxLength, {1, 2, 3}, {4, 5, 6}, bottomLeftCorner);

  _configuration->cubeGaussObjects.push_back(cubeGauss);
  _configuration->calcSimulationBox();

  std::array<double, 3> expectedBoxMin;
  std::array<double, 3> expectedBoxMax;
  for (unsigned int i = 0; i < expectedBoxMax.size(); ++i) {
    expectedBoxMax[i] = bottomLeftCorner[i] + boxLength[i];
    // box min is at least 0 (see initialization of _config)
    expectedBoxMin[i] = std::min(bottomLeftCorner[i], 0.);
  }
  EXPECT_THAT(_configuration->boxMin.value, ElementsAreArray(expectedBoxMin));
  EXPECT_THAT(_configuration->boxMax.value, ElementsAreArray(expectedBoxMax));
}

TEST_F(MDFlexConfigTest, UniformBoxMinMax) {
  size_t numParticles = 10;
  std::array<double, 3> boxLength{4, 5, 6};
  std::array<double, 3> bottomLeftCorner{-1, 0, 1};
  CubeUniform cubeUniform(zero, 0, 1, 1, 1, numParticles, boxLength, bottomLeftCorner);

  _configuration->cubeUniformObjects.push_back(cubeUniform);
  _configuration->calcSimulationBox();

  std::array<double, 3> expectedBoxMin;
  std::array<double, 3> expectedBoxMax;
  for (unsigned int i = 0; i < expectedBoxMax.size(); ++i) {
    expectedBoxMax[i] = bottomLeftCorner[i] + boxLength[i];
    // box min is at least 0 (see initialization of _config)
    expectedBoxMin[i] = std::min(bottomLeftCorner[i], 0.);
  }
  EXPECT_THAT(_configuration->boxMin.value, ElementsAreArray(expectedBoxMin));
  EXPECT_THAT(_configuration->boxMax.value, ElementsAreArray(expectedBoxMax));
}
