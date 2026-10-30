/**
 * @file TuningManagerTest.h
 * @author muehlhaeusser
 * @date 18.04.2026
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/tuning/Configuration.h"

class TuningManagerTest : public AutoPasTestBase {
 public:
  TuningManagerTest() = default;
  ~TuningManagerTest() override = default;

  /**
   * Tests that ending a tuning phase with a rejected configuration is handled correctly.
   *
   * Mimics a TuningManager with one AutoTuner trialling four configurations. The second and the final configurations
   * will be rejected.
   *
   * @param rejectIndefinitely Run the test, passing this rejectIndefinitely value to rejectConfiguration
   */
  void testEndingTuningPhaseWithRejectedConfig(bool rejectIndefinitely) const;
};
