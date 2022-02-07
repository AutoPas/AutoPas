/**
 * @file MixedBoundaryConditionTest.h
 * @author S. J. Newcome
 * @date 21/01/2022
 */
#pragma once

#include <gtest/gtest.h>
//#include <gmock/gmock-matchers.h>

#include "AutoPasTestBase.h"

/**
 * Test class for the addition of mixed boundary conditions to RegularGridDecomposition
 */
class MixedBoundaryConditionTest : public AutoPasTestBase {
 public:
  /**
   * Constructor.
   */
  MixedBoundaryConditionTest() = default;
};

/**
 * Parameterized test case for reflective boundary conditions in RegularGridDecomposition
 */
class ReflectiveBoundaryConditionTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple</*position*/ std::array<double, 3>,
                                                      /*velocity*/ std::array<double, 3>,
                                                      /*isReflected*/ std::array<bool, 3>>> {
 public:
  /**
   * Constructor.
   */
  ReflectiveBoundaryConditionTest() = default;
};