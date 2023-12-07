/**
 * @file NumThreadGuard.h
 * @author C. Menges
 * @date 14.04.2019
 */

#pragma once

#include <autopas/utils/WrapOpenMP.h>

/**
 * NumThreadGuard sets current number of threads to newNum and resets number of threads during destruction.
 */
class NumThreadGuard final {
 public:
  /**
   * Construct a new NumThreadGuard object and sets current number of threads to newNum.
   * @param newNum new number of threads
   */
  explicit NumThreadGuard(const int newNum) : numThreadsBefore(autopas::autopas_get_max_threads()) {
    autopas::autopas_set_num_threads(newNum);
  }

  /**
   * Destroy the NumThreadGuard object and reset number of threads.
   */
  ~NumThreadGuard() { autopas::autopas_set_num_threads(numThreadsBefore); }

  /**
   * delete copy constructor.
   */
  NumThreadGuard(const NumThreadGuard &) = delete;

  /**
   * delete copy assignment constructor.
   * @return deleted, so not important
   */
  NumThreadGuard &operator=(const NumThreadGuard &) = delete;

 private:
  int numThreadsBefore;
};
