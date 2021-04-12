/**
 * @file ArcherShouldReportThis.cpp
 * @author seckler
 * @date 25.05.18
 */
#ifdef AUTOPAS_OPENMP

#include <gtest/gtest.h>
#include <omp.h>

#include <atomic>
#include <numeric>

class Container {
  enum class ValidityState : unsigned char {
    invalid = 0,
    valid = 1,
  };

  std::vector<int> vec{};
  std::atomic<ValidityState> _isValid{ValidityState::invalid};

  auto generate() {
    vec = {1, 2, 3};
    _isValid = ValidityState::valid;
  }

 public:
  auto begin_critical() {
#pragma omp critical
    if (_isValid == ValidityState::invalid) {
      generate();
    }

    return vec.begin() + omp_get_thread_num();
  }

  auto begin_single() {
#pragma omp single
    if (_isValid == ValidityState::invalid) {
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
  for (auto iter = c.begin_single(); iter < c.end(); iter += omp_get_num_threads()) {
    *iter = *iter * 2;
    res += *iter;
  }

  int reference_res{0};
  for (auto iter = c.begin_single(); iter < c.end(); iter += 1) {
    reference_res += *iter;
  }
  EXPECT_GT(res, 0);
  EXPECT_EQ(res, reference_res);
}

TEST(ArcherShouldNotReportThis, testGenerateUsingCriticalInBegin) {
  Container c;
  int res{0};
#pragma omp parallel reduction(+ : res)
  for (auto iter = c.begin_critical(); iter < c.end(); iter += omp_get_num_threads()) {
    *iter = *iter * 2;
    res += *iter;
  }

  int reference_res{0};
  for (auto iter = c.begin_critical(); iter < c.end(); iter += 1) {
    reference_res += *iter;
  }
  EXPECT_GT(res, 0);
  EXPECT_EQ(res, reference_res);
}

#endif