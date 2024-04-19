/**
 * @file WrapOpenMP.h
 * @author F. Gratl
 * @date 4/20/18
 *
 * @details
 * Provide non-OpenMP versions of the most common OpenMP function calls,
 * so that they don't have to be wrapped in ifdef-s every time.
 *
 * Proper wrapper and renaming necessary, because of -fopenmp-simd handling of
 * gcc.
 *
 * Extend when necessary.
 */

#pragma once

#if defined(AUTOPAS_USE_OPENMP)
#include <omp.h>

#include <cstddef>  // for size_t
#include <vector>

#else
#include "ExceptionHandler.h"
#endif

namespace autopas {

#if defined(AUTOPAS_USE_OPENMP)

/**
 * Helper macro to stringify arguments and use them as a pragma
 * This is necessary to combine literals and arguments to a string argument for _Pragma()
 */
#define AUTOPAS_DO_PRAGMA(x) _Pragma(#x)

/**
 * Wrapper macro to replace "#pragma omp" that can be (de)activated.
 * @param args Anything one can pass to `#pragma omp`
 *
 * @note Only one argument is supported. OpenMP clauses must thus be separated by a space.
 */
#define AUTOPAS_OPENMP(args) AUTOPAS_DO_PRAGMA(omp args)

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
 * Wrapper for omp_set_num_threads().
 * @param n New max number of threads.
 */
inline void autopas_set_num_threads(int n) { omp_set_num_threads(n); }

/**
 * Wrapper for omp_set_schedule().
 * Sets the scheduling kind and chunk size used by schedule(runtime).
 * @param kind the scheduling kind to use
 * @param chunkSize the chunk size to use
 */
inline void autopas_set_schedule(omp_sched_t kind, int chunkSize) { omp_set_schedule(kind, chunkSize); }

/**
 * Wrapper for omp_get_schedule().
 * Puts the values of OpenMP's scheduling runtime variables at the given pointers.
 */
inline void autopas_get_schedule(omp_sched_t *kind, int *chunkSize) { omp_get_schedule(kind, chunkSize); }

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
  AutoPasLock(AutoPasLock &&) noexcept { omp_init_lock(&_lock); }

  /**
   * Copy constructor
   */
  AutoPasLock(const AutoPasLock &) { omp_init_lock(&_lock); }

  /**
   * Assignment operator
   * @return reference to this object after copy
   */
  AutoPasLock &operator=(AutoPasLock) = delete;

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
#pragma omp declare reduction( \
        vecMerge : std::vector<size_t> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction( \
        vecMerge : std::vector<double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

#else

/**
 * Empty macro to throw away any arguments.
 */
#define AUTOPAS_OPENMP(args)

/**
 * Wrapper for omp_sched_t, same as in OpenMP's omp.h
 */
typedef enum omp_sched_t {
  omp_sched_static = 1,
  omp_sched_dynamic = 2,
  omp_sched_guided = 3,
  omp_sched_auto = 4,
  omp_sched_monotonic = 0x80000000U
} omp_sched_t;

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
 * Wrapper for omp_set_num_threads().
 * Does nothing when OpenMP is disabled.
 */
inline void autopas_set_num_threads(int /* n */) {}

/**
 * Wrapper for omp_set_schedule().
 * @param kind the scheduling kind to use
 * @param chunkSize the chunk size to use
 */
inline void autopas_set_schedule(OpenMPKindOption /* kind */, int /* chunkSize */) {}

/**
 * Wrapper for omp_get_schedule().
 * Puts the values of OpenMP's scheduling runtime variables at the given pointers.
 */
inline void autopas_get_schedule(omp_sched_t *kind, int *chunkSize) {}

/**
 * AutoPasLock for the sequential case.
 */
class AutoPasLock {
 public:
  /**
   * Default constructor
   */
  AutoPasLock() { _locked = false; }

  /**
   * Move Constructor
   */
  AutoPasLock(AutoPasLock &&) noexcept { _locked = false; }

  /**
   * Copy constructor
   */
  AutoPasLock(AutoPasLock &) { _locked = false; }

  /**
   * Assignment operator
   * @return reference to this object after copy
   */
  AutoPasLock &operator=(AutoPasLock) = delete;

  /**
   * Destructor
   */
  ~AutoPasLock() {
    if (_locked) {
      utils::ExceptionHandler::exception("AutoPasLocked destroyed in locked state.");
    }
  }

  /**
   * Acquire the lock.
   */
  void lock() {
    if (_locked) {
      utils::ExceptionHandler::exception("Tried to acquire a locked lock.");
    }
    _locked = true;
  }

  /**
   * Release the lock.
   */
  void unlock() {
    if (not _locked) {
      utils::ExceptionHandler::exception("Tried to release an unlocked lock.");
    }
    _locked = false;
  }

 private:
  // true if locked, false if unlocked
  bool _locked;
};

#endif

// These properties are needed because we use AutoPasLock in vectors on which we call resize().
static_assert(std::is_default_constructible_v<AutoPasLock>, "AutoPasLock needs to be default constructible!");
static_assert(std::is_move_constructible_v<AutoPasLock>, "AutoPasLock needs to be move constructible!");

}  // namespace autopas