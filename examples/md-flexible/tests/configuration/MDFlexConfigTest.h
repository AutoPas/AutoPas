/**
 * @file MDFlexConfigTest.h
 * @author F. Gratl
 * @date 04/06/2020
 */

#pragma once
#include <gtest/gtest.h>

#include <array>

#include "AutoPasTestBase.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/configuration/objects/CubeGauss.h"
#include "src/configuration/objects/CubeGrid.h"
#include "src/configuration/objects/CubeUniform.h"
#include "src/configuration/objects/Sphere.h"

class MDFlexConfigTest : public AutoPasTestBase {
 public:
  MDFlexConfigTest() : AutoPasTestBase() { }

 protected:
  std::shared_ptr<MDFlexConfig> _configuration;
  double _interactionLength;

  /**
   * Utility zero-vector.
   */
  std::array<double, 3> zero{0, 0, 0};
};
