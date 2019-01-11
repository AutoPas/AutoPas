/**
 * @file WrapOpenMP.h
 * @author F. Gratl
 * @date 4/20/18
 */

#pragma once

/**
 * Provide non-OpenMP versions of the most common OpenMP function calls,
 * so that they don't have to be wrapped in ifdef-s every time.
 *
 * Proper wrapper and renaming necessary, because of -fopenmp-simd handling of
 * gcc.
 *
 * Extend when necessary.
 */

#include <cassert>
#if defined(AUTOPAS_OPENMP)
#include <omp.h>
#endif

namespace autopas {

#if defined(AUTOPAS_OPENMP)

/**
 * Wrapper for omp_get_thread_num().
 * @return Id of the current thread.
 */
inline int autopas_get_thread_num() { return omp_get_thread_num(); }

/**
 * Wrapper for omp_get_num_thread().
 * @return Number of currently active threads.
 */
inline int autopas_get_num_threads() { return omp_get_num_threads(); }

/**
 * Wrapper for omp_get_max_threads().
 * @return Number of threads that can be activated.
 */
inline int autopas_get_max_threads() { return omp_get_max_threads(); }

/**
 * AutoPasLock for the openmp case, this wraps a omp_lock_t object. To make it copyable, etc.
 */
class AutoPasLock {
 public:
  /**
   * Default constructor
   */
  AutoPasLock() { omp_init_lock(&_lock); }

  /**
   * Move Constructor
   */
  AutoPasLock(AutoPasLock&&) noexcept {
    omp_init_lock(&_lock);
  }

  /**
   * Copy constructor
   */
  AutoPasLock(const AutoPasLock&) {
    omp_init_lock(&_lock);
  }

  /**
   * Assignment operator
   * @return reference to this object after copy
   */
  AutoPasLock& operator=(AutoPasLock) = delete;

  /**
   * Destructor
   */
  ~AutoPasLock() { omp_destroy_lock(&_lock); }

  /**
   * Acquire the lock.
   */
  void lock() { omp_set_lock(&_lock); }

  /**
   * Release the lock.
   */
  void unlock() { omp_unset_lock(&_lock); }

 private:
  omp_lock_t _lock;
};

/**
 * Custom reductions:
 */
// reduction for merging vectors: {1,2} + {2,3} -> {1,2,2,3}
#pragma omp declare reduction(vecMerge : std::vector<size_t> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(vecMerge : std::vector<double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

#else

/**
 * Dummy for omp_set_lock() when no OpenMP is available.
 * @return Always 0.
 */
inline int autopas_get_thread_num() { return 0; }

/**
 * Dummy for omp_get_num_threads() when no OpenMP is available.
 * @return Always 1.
 */
inline int autopas_get_num_threads() { return 1; }

/**
 * Dummy for omp_get_max_threads() when no OpenMP is available.
 * @return Always 1.
 */
inline int autopas_get_max_threads() { return 1; }

/**
 * AutoPasLock for the sequential case, that uses an enum for checking the base functionality.
 */
class AutoPasLock {
 public:
  /**
   * Default constructor
   */
  AutoPasLock() { _lock = unlocked; }

  /**
   * Move Constructor
   */
  AutoPasLock(AutoPasLock&&) noexcept {
    _lock = unlocked;
  }

  /**
   * Copy constructor
   */
  AutoPasLock(AutoPasLock&) {
    _lock = unlocked;
  }

  /**
   * Assignment operator
   * @return reference to this object after copy
   */
  AutoPasLock& operator=(AutoPasLock) = delete;

  /**
   * Destructor
   */
  ~AutoPasLock() {
    assert(_lock == unlocked);
  }

  /**
   * Acquire the lock.
   */
  void lock() {
    assert(_lock == unlocked);
    _lock = locked;
  }

  /**
   * Release the lock.
   */
  void unlock() {
    assert(_lock == locked);
    _lock = unlocked;
  }

 private:
  // lock: 0 means unlocked, 1 locked.
  enum {unlocked, locked} _lock;
};

#endif

}  // namespace autopas