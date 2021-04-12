/**
 * @file ArcherShouldReportThis.cpp
 * @author seckler
 * @date 25.05.18
 */
#ifdef AUTOPAS_OPENMP

#include <gtest/gtest.h>
#include <omp.h>

#include <numeric>

class Container {
  std::vector<int> vec{0, 0, 0};
  bool generated = false;

  auto generate() {
    vec = {1, 2, 3};
    std::iota(vec.begin(), vec.end(), 1);
    generated = true;
  }

 public:
  auto begin() {
#pragma omp single
    if (not generated) {
      generate();
    }

    return vec.begin() + omp_get_thread_num();
  }

  auto end() { return vec.end(); }
};

TEST(ArcherShouldNotReportThis, testGenerateUsingSingleInBegin) {
  Container c;
  int res{0};
#pragma omp parallel reduction(+ : res)
  for (auto iter = c.begin(); iter < c.end(); iter += omp_get_num_threads()) {
    *iter = *iter * 2;
    res += *iter;
  }

  int reference_res{0};
  for (auto iter = c.begin(); iter < c.end(); iter += 1) {
    reference_res += *iter;
  }
  EXPECT_GT(res, 0);
  EXPECT_EQ(res, reference_res);
}

#endif