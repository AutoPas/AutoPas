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

  /**
   * Derive expected positions and velocities according to the setup.
   * @param particlePositions
   * @param particleVelocities
   * @param boxMin global box min
   * @param boxMax global box max
   * @param sigma sigma of particle used
   * @param interactionLength
   * @param boundaryConditions
   * @return
   */
  static auto setUpExpectations(const std::vector<std::array<double, 3>> &particlePositions,
                                const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double sigma,
                                double interactionLength,
                                const std::array<options::BoundaryTypeOption, 3> &boundaryConditions);

  void testFunction(const std::vector<std::array<double, 3>> &particlePositions,
                    const std::array<options::BoundaryTypeOption, 3> &boundaryConditions);
};
