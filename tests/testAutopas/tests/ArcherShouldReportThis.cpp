/**
 * @file ArcherShouldReportThis.cpp
 * @author seckler
 * @date 25.05.18
 */

#include <gtest/gtest.h>

namespace ArcherShouldReportThis {

constexpr size_t N = 100;

// this test is here just to test if archer is working.
// and is thus currently disabled.
// remove the "DISABLED_" part to enable it again.
#ifdef AUTOPAS_OPENMP
TEST(DISABLED_ArcherShouldReportThis, test) {
  int a[N];
  for (int i = 0; i < N; i++) {
    a[i] = i % 7;
  }

#pragma omp parallel for
  for (int i = 0; i < N - 1; i++) {
    a[i] = a[i + 1];
  }
  std::cout << "a[80]" << a[80] << std::endl;
}
#endif
} // end namespace ArcherShouldReportThis
