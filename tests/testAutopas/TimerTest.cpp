/**
 * @file TimerTest.cpp
 * @author seckler
 * @date 17.04.18
 */

#include <gtest/gtest.h>
#include "autopas/utils/Timer.h"

TEST(TimerTest, testTimer) {
  autopas::utils::Timer timer;
  timer.start();
  auto time = timer.stop();
  ASSERT_LE(time, 0.01);

  timer.start();
  // sleeping for 100ms
  usleep(5e4);
  time = timer.stop();
  // time should be close to 10ms
  ASSERT_NEAR(time, .05, .03);
}