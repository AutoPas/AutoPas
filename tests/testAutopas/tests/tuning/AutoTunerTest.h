/**
 * @file AutoTunerTest.h
 * @author F. Gratl
 * @date 8/10/18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/tuning/Configuration.h"

class AutoTunerTest : public AutoPasTestBase {
 public:
  AutoTunerTest() = default;
  ~AutoTunerTest() override = default;

  /**
   * Tests that ending a tuning phase with a reject configuration is handled correctly.
   *
   * Mimics an AutoTuner trialling four configurations. The second and the final configurations will be rejected. After
   * the tuning phase is completed, the test should get the correct configuration upon calling `getNextConfig`
   *
   * @param rejectIndefinitely Run the test, parsing this rejectIndefinitely value to rejectConfig
   */
  void testEndingTuningPhaseWithRejectedConfig(bool rejectIndefinitely) const;
};
