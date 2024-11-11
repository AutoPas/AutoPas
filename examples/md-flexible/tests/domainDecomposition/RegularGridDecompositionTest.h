/**
 * @file RegularGridDecompositionTest.h
 * @author J. Körner
 * @date 27/05/21
 */
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "src/options/BoundaryTypeOption.h"

/**
 * Test class for the RegularGridDecomposition domain decomposition class.
 */
class RegularGridDecompositionTest : public AutoPasTestBase,
                                     public ::testing::WithParamInterface<options::BoundaryTypeOption> {
 public:
  /**
   * Constructor.
   */
  RegularGridDecompositionTest() = default;
};
