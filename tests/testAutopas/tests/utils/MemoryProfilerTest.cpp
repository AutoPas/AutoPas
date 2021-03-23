/**
 * @file MemoryProfilerTest.cpp.cc
 * @author F. Gratl
 * @date 23.03.21
 */

#include "MemoryProfilerTest.h"

#include "autopas/utils/MemoryProfiler.h"

/**
 * Make sure the memory profiler finds any memory usage.
 */
TEST(MemoryProfilerTest, currentMemoryUsage) {
  EXPECT_GT(autopas::memoryProfiler::currentMemoryUsage(), 0)
      << "Memory profiler did not find any memory usage! Typically this test should report ~3300kb.";
}