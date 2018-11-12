/**
 * @file AutoTunerTest.h
 * @author F. Gratl
 * @date 8/10/18
 */

#pragma once

#include <gtest/gtest.h>
#include <autopas/selectors/AutoTuner.h>
#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

class AutoTunerTest : public AutoPasTestBase {
 public:
  AutoTunerTest() = default;
  ~AutoTunerTest() = default;

  // @todo when SoA and Newton 3 are tuneable extend this.
  void testTune(autopas::DataLayoutOption dataLayoutOption);
};
