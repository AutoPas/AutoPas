/**
 * @file TraceTimerTest.cpp
 * @author muehlhaeusser
 * @date 12.02.2026
 */

#include "TraceTimerTest.h"

#include <gtest/gtest.h>

#include <chrono>
#include <thread>

TEST_F(TraceTimerTest, TestStartStop) {
  _timer.start();

  // Sleep for a tiny amount to ensure non-zero time passes
  std::this_thread::sleep_for(std::chrono::milliseconds(10));

  long elapsed = _timer.stop();

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  // We expect the timer to actually measure time.
  EXPECT_GT(elapsed, 0) << "Timer should return > 0 ns when Trace is active.";
  EXPECT_GT(_timer.getTotalTime(), 0) << "Total time should be updated.";

#else
  // We expect the timer to be a 'no-op' and return 0.
  EXPECT_EQ(elapsed, 0) << "Timer should return 0 when Trace is inactive.";
  EXPECT_EQ(_timer.getTotalTime(), 0) << "Total time should remain 0.";
#endif
}

TEST_F(TraceTimerTest, TestReset) {
  _timer.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(1));
  _timer.stop();

  // Verify we have time (if active)
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  ASSERT_GT(_timer.getTotalTime(), 0);
#endif

  _timer.reset();

  // Verify reset worked
  EXPECT_EQ(_timer.getTotalTime(), 0);
}