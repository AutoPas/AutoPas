/**
 * @file TimerTest.cpp
 * @author seckler
 * @date 17.04.18
 */

#include <gtest/gtest.h>

#include <chrono>
#include <thread>

#include "autopas/utils/Timer.h"

namespace TimerTest {

TEST(DISABLED_TimerTest, testTimer) {
  autopas::utils::Timer timer;
  timer.start();
  auto time = timer.stop();
  ASSERT_LE(time, 0.01);

  timer.start();
  // sleeping for 0.05s
  using namespace std::chrono_literals;
  std::this_thread::sleep_for(50ms);

  time = timer.stop();
  // time should be close to 50ms
  ASSERT_NEAR(time, .05, .03);
}
} // end namespace TimerTest
