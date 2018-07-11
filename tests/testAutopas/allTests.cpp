/**
 * @file allTests.cpp
 * @author seckler
 * @date 18.01.18
 */

#include <gtest/gtest.h>

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
