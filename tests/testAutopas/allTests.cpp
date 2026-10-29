/**
 * @file allTests.cpp
 * @author seckler
 * @date 18.01.18
 */

#include <gtest/gtest.h>

#ifdef AUTOPAS_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  // damit nur meine Tests durchgelaufen werden
  const auto result = RUN_ALL_TESTS();

#ifdef AUTOPAS_ENABLE_KOKKOS
  // Kokkos is initialized lazily by KokkosTestBase, so only finalize it if a Kokkos test actually ran.
  if (Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
#endif

  return result;
}
