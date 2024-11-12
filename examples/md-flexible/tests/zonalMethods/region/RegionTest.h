
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "src/TypeDefinitions.h"

/**
 * Test class for the Zone classes
 */
class RegionTest : public AutoPasTestBase {
 public:
  using AutoPasType = autopas::AutoPas<ParticleType>;

  /**
   * Constructor.
   */
  RegionTest() : AutoPasTestBase() {}

 protected:
  // initialize autopas container
  void initContainer(AutoPasType &autopas, std::vector<ParticleType> particles);

  AutoPasType _autopas;
};
