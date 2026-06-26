/**
 * @file TimerTest.cpp
 * @author S. Newcome
 * @date 26/06/2026
 */

#include "TimerTest.h"

#include <gtest/gtest.h>

#include <cctype>
#include <chrono>
#include <string>
#include <thread>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Timer.h"

using autopas::utils::ExceptionBehavior;
using autopas::utils::ExceptionHandler;
using autopas::utils::Timer;

void TimerTest::SetUp() { ExceptionHandler::setBehavior(ExceptionBehavior::throwException); }

/**
 * A freshly constructed Timer has neither accumulated total time nor lap time.
 */
TEST_F(TimerTest, Constructor) {
  Timer timer;
  EXPECT_EQ(timer.getTotalTime(), 0);
  EXPECT_EQ(timer.getLapTime(), 0);
}

/**
 * The normal usage path (start() followed by stop()) on a fresh timer must not throw.
 */
TEST_F(TimerTest, StartThenStopDoesNotThrow) {
  Timer timer;
  timer.start();
  EXPECT_NO_THROW(timer.stop());
}

/**
 * Calling start() on a timer that is already running is a misuse and must throw.
 */
TEST_F(TimerTest, StartWhileRunningThrows) {
  Timer timer;
  timer.start();
  EXPECT_THROW(timer.start(), ExceptionHandler::AutoPasException);
}

/**
 * Calling stop() on a timer that was never started is a misuse and must throw.
 */
TEST_F(TimerTest, StopWithoutStartThrows) {
  Timer timer;
  EXPECT_THROW(timer.stop(), ExceptionHandler::AutoPasException);
}

/**
 * stop() returns the elapsed time in nanoseconds since the last start(). As sleep_for only provides
 * a lower bound on the slept duration, the returned value is only checked to be at least kSleep.
 */
TEST_F(TimerTest, StopReturnsElapsedLowerBound) {
  Timer timer;
  timer.start();
  std::this_thread::sleep_for(kSleep);
  const long elapsed = timer.stop();
  EXPECT_GE(elapsed, kSleepNanos);
}

/**
 * stop() accumulates the measured duration into both the total time and the lap time. Both must be
 * at least the slept duration (lower bound only).
 */
TEST_F(TimerTest, StopAccumulatesTotalAndLap) {
  Timer timer;
  timer.start();
  std::this_thread::sleep_for(kSleep);
  timer.stop();
  EXPECT_GE(timer.getTotalTime(), kSleepNanos);
  EXPECT_GE(timer.getLapTime(), kSleepNanos);
}

/**
 * reset() clears both the total time and the lap time. Tested with addTime() to keep the assertion
 * independent of clock timing.
 */
TEST_F(TimerTest, Reset) {
  Timer timer;
  timer.addTime(1234);
  timer.reset();
  EXPECT_EQ(timer.getTotalTime(), 0);
  EXPECT_EQ(timer.getLapTime(), 0);
}

/**
 * resetLap() clears only the lap time and leaves the accumulated total time untouched. Tested with
 * addTime() to keep the assertion independent of clock timing.
 */
TEST_F(TimerTest, ResetLap) {
  Timer timer;
  timer.addTime(1234);
  timer.resetLap();
  EXPECT_EQ(timer.getLapTime(), 0);
  EXPECT_EQ(timer.getTotalTime(), 1234);
}

/**
 * addTime() adds the given amount of nanoseconds to both the total and lap time, and repeated calls
 * accumulate. This is fully deterministic as it does not depend on the clock.
 */
TEST_F(TimerTest, AddTime) {
  Timer timer;

  timer.addTime(100);
  EXPECT_EQ(timer.getTotalTime(), 100);
  EXPECT_EQ(timer.getLapTime(), 100);

  timer.addTime(50);
  EXPECT_EQ(timer.getTotalTime(), 150);
  EXPECT_EQ(timer.getLapTime(), 150);
}