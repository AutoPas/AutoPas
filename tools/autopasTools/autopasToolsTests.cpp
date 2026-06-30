/**
 * @file autopasToolsTests.cpp
 * @author hmeyran
 * @date 30.06.2026
 */

#include <gtest/gtest.h>

/**
 * GTest entry point for autopasTools tests.
 */
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
