/**
 * @file ArcherShouldReportThis.cpp
 * @author seckler
 * @date 25.05.18
 */

#include <gtest/gtest.h>

#define N 100

#ifdef AUTOPAS_OPENMP
TEST(ArcherShouldReportThis, test) {
  int a[N];
  for(int i=0;i<N;i++){
    a[i]=i%7;
  }

#pragma omp parallel for
  for (int i = 0; i < N - 1; i++) {
    a[i] = a[i + 1];
  }
  std::cout << "a[80]" << a[80] << std::endl;
}
#endif