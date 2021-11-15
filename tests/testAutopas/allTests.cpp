/**
 * @file allTests.cpp
 * @author seckler
 * @date 18.01.18
 */

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  // damit nur meine Tests durchgelaufen werden
  Kokkos::initialize(argc, argv);
  auto result = RUN_ALL_TESTS();
  Kokkos::finalize();

  return result;
}
