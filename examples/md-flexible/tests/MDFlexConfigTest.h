/**
 * @file MDFlexConfigTest.h
 * @author F. Gratl
 * @date 04/06/2020
 */

#pragma once
#include <gtest/gtest.h>
#include <src/parsing/MDFlexConfig.h>

#include <array>

#include "AutoPasTestBase.h"
#include "src/Objects/CubeGauss.h"
#include "src/Objects/CubeGrid.h"
#include "src/Objects/CubeUniform.h"
#include "src/Objects/Sphere.h"

class MDFlexConfigTest : public AutoPasTestBase {
 public:
  MDFlexConfigTest() : AutoPasTestBase() {
    _config.cutoff.value = 1;
    _config.verletSkinRadius.value = .5;
    _config.boxMin.value = {0., 0., 0.};
    _config.boxMax.value = {1., 1., 1.};

    _interactionLength = _config.cutoff.value + _config.verletSkinRadius.value;
  }

 protected:
  MDFlexConfig _config;
  double _interactionLength;

  /**
   * Utility zero-vector.
   */
  std::array<double, 3> zero{0, 0, 0};
};