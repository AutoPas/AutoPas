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
  MDFlexConfigTest() : AutoPasTestBase() {
    // std::string arguments =
    //"md-flexible --yaml-filename " + std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";
    //
    //_configuration = std::make_shared<MDFlexConfig>(3, reinterpret_cast<char **>(arguments.data()));
    //
    //_configuration->cutoff.value = 1;
    //_configuration->verletSkinRadius.value = .5;
    //_configuration->boxMin.value = {0., 0., 0.};
    //_configuration->boxMax.value = {1., 1., 1.};
    //
    //_interactionLength = _configuration->cutoff.value + _configuration->verletSkinRadius.value;
  }

 protected:
  std::shared_ptr<MDFlexConfig> _configuration;
  double _interactionLength;

  /**
   * Utility zero-vector.
   */
  std::array<double, 3> zero{0, 0, 0};
};
