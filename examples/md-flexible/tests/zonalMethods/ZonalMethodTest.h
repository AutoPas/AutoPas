
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "src/zonalMethods/ZonalMethod.h"

/**
 * Test class for the Zone classes
 */
class ZonalMethodTest : public ZonalMethod, public AutoPasTestBase {
 public:
  /**
   * Constructor.
   */
  ZonalMethodTest() : ZonalMethod(1), AutoPasTestBase() {}

  void collectParticles(AutoPasType &autoPasContainer) override;
};
