/**
 * @file TimerTest.h
 * @author S. Newcome
 * @date 26/06/2026
 */

#pragma once

#include <chrono>

#include "AutoPasTestBase.h"

/**
 * Test fixture for the Timer class.
 *
 * It forces the global ExceptionHandler behavior to throwException in SetUp so that the tests
 * verifying misuse of start()/stop() can rely on exceptions being thrown, independent of the order
 * in which gtest runs the test cases (other test suites mutate this global state).
 */
class TimerTest : public AutoPasTestBase {
 protected:
  void SetUp() override;

  /**
   * Duration the timing tests sleep for before stopping the timer.
   *
   * std::this_thread::sleep_for only guarantees a *lower* bound: the calling thread is blocked for
   * *at least* this duration, but possibly (much) longer under load, which is common in CI. Tests
   * therefore only ever assert that the measured time is greater than or equal to this value and
   * never assert a tight upper bound.
   */
  static constexpr auto kSleep = std::chrono::milliseconds(20);

  /**
   * kSleep expressed in nanoseconds, matching the unit returned by Timer::stop() / getTotalTime().
   */
  static constexpr long kSleepNanos = std::chrono::duration_cast<std::chrono::nanoseconds>(kSleep).count();
};
