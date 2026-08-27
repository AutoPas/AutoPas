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
   * Tests if thread count tuning selected a valid option out of a given set and (optionally) the the expected number of
   * threads (e.g. minimum for a small simulation size, maximum for a large simulation size, or max threads if 0).
   *
   * @param boxMax The size of the simulation domain along each side
   * @param threadCountOptions The possible thread counts in the configuration
   * @param expectedSelectedThreadCount Optional expected number of threads of the configuration after tuning has
   * completed (0 = max threads)
   */
  void testThreadCountTuningWithBoxMax(const size_t boxMax, const std::set<int> &threadCountOptions,
                                       int expectedSelectedThreadCount = -1) const;
};
