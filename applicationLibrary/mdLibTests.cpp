/**
 * @file mdLibTests.cpp
 * @author S. Newcome
 * @date 12/05/2023
 */

#include <gtest/gtest.h>

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  int output = RUN_ALL_TESTS();
  return output;
}
