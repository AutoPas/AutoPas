/**
 * @file NumThreadGuard.h
 * @author C. Menges
 * @date 14.04.2019
 */

#pragma once

#include "autopas/utils/WrapOpenMP.h"

/**
 * Thread setting to guard with \ref NumThreadGuard
 */
enum RestoreType {
  /**
   * Guard \ref autopas::autopas_set_num_threads
   */
  USED_THREADS,
  /**
   * Guard \ref autopas::autopas_set_tuned_num_threads
   */
  TUNED_THREADS,
};

/**
 * NumThreadGuard sets current number of threads to newNum and resets number of threads during destruction.
 */
class NumThreadGuard final {
 public:
  /**
   * Construct a new NumThreadGuard object and sets current number of threads to newNum.
   * @param newNum new number of threads
   * @param restoreType Type of thread setting to set/restore (Defaults to USED_THREADS)
   */
  explicit NumThreadGuard(const int newNum, const RestoreType restoreType = USED_THREADS) : _restoreType(restoreType) {
    if (_restoreType == USED_THREADS) numThreadsBefore = autopas::autopas_get_num_threads();
    if (_restoreType == TUNED_THREADS)
      numThreadsBefore = autopas::autopas_get_tuned_num_threads();
    else
      std::runtime_error("Unknown RestoreType!");
    set_num_threads(newNum);
  }

  /**
   * Destroy the NumThreadGuard object and reset number of threads.
   */
  ~NumThreadGuard() { set_num_threads(numThreadsBefore); }

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
  void set_num_threads(const int newNum) {
    if (_restoreType == USED_THREADS) autopas::autopas_set_num_threads(newNum);
    if (_restoreType == TUNED_THREADS)
      autopas::autopas_set_tuned_num_threads(newNum);
    else
      std::runtime_error("Unknown RestoreType!");
  }

  int numThreadsBefore;
  RestoreType _restoreType;
};
