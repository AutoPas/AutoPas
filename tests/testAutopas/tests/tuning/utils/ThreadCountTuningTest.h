/**
 * @file ThreadTuning.h
 * @author R. Horn
 * @date 06/03/2026
 */

#pragma once

#include "AutoPasTestBase.h"

class ThreadCountTuningTest : public AutoPasTestBase {
 public:
  /**
   * Tests if thread count tuning selected the expected number of threads for a given simulation size
   * 
   * @param boxMax The size of the simulation domain along each side
   * @param threadCountOptions The possible thread counts in the configuration
   * @param expectedSelectedThreadCount The expected number of threads of the configuration after tuning has completed
   */
  void testThreadCountTuningWithBoxMax(const size_t boxMax, const std::set<int> &threadCountOptions, const size_t expectedSelectedThreadCount) const;
};
