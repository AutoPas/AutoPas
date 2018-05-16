/**
 * @file TimerTest.cpp
 * @author seckler
 * @date 17.04.18
 */

#include <gtest/gtest.h>
#include "utils/Timer.h"

TEST(TimerTest, testTimer) {
  autopas::utils::Timer timer;
  timer.start();
  auto time = timer.stop();
  ASSERT_LE(time, 0.01);

  timer.start();
  // sleeping for 100ms
  usleep(1e5);
  time = timer.stop();
  // time should be close to 10ms
  ASSERT_NEAR(time, 1e-5, 5e-4);
}