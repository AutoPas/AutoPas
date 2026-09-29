/**
 * @file GeneratorsTest.h
 * @author N. Fottner
 * @date 02/08/19
 */
#pragma once
#include <gtest/gtest.h>

#include <array>

#include "AutoPasTestBase.h"
#include "src/configuration/objects/CubeClosestPacked.h"

class GeneratorsTest : public AutoPasTestBase {
 public:
  GeneratorsTest() = default;

 protected:
  double epsilon{1.0};
  double sigma{1.0};
  double cutoff{1.};
  std::array<double, 3> boxmin{{0., 0., 0.}};
  std::array<double, 3> boxmax{{5., 5., 5.}};

  /**
   * Test scenario structure for CubeGrid constructed with target density.
   */
  struct GridDensityScenario {
    std::array<size_t, 3> particlesPerDim;
    std::array<double, 3> bottomLeft;
    double density;
    bool centered;
    std::string description;
  };

  /**
   * Test scenario structure for CubeUniform constructed with target density.
   */
  struct UniformDensityScenario {
    std::array<double, 3> boxLength;
    std::array<double, 3> bottomLeft;
    double density;
    std::string description;
  };

  /**
   * Test scenario structure for CubeClosestPacked constructed with target density.
   */
  struct ClosestPackedDensityScenario {
    CubeClosestPacked::LatticeStructure structure;
    std::array<double, 3> boxLength;
    std::array<double, 3> bottomLeft;
    double density;
    bool centered;
    std::string description;
  };
};