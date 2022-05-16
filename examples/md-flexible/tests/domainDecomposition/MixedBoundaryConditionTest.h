/**
 * @file MixedBoundaryConditionTest.h
 * @author S. J. Newcome
 * @date 21/01/2022
 */
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "src/options/BoundaryTypeOption.h"

/**
 * Test class for the addition of mixed boundary conditions to RegularGridDecomposition
 */
class MixedBoundaryConditionTest : public AutoPasTestBase {
 public:
  /**
   * Constructor.
   */
  MixedBoundaryConditionTest() = default;

  static auto setUpExpectations(const std::vector<std::array<double, 3>> &particlePositions,
                                const std::vector<std::array<double, 3>> &particleVelocities,
                                const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                const double reflectionSkin, const double interactionLength,
                                const std::array<options::BoundaryTypeOption, 3> &boundaryConditions);

  void testFunction(const std::vector<std::array<double, 3>> &particlePositions,
                    const std::vector<std::array<double, 3>> &particleVelocities,
                    const std::array<options::BoundaryTypeOption, 3> &boundaryConditions);
};
