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
 * Wrapper for omp_lock_t.
 */
typedef omp_lock_t autopas_lock_t;

/**
 * Wrapper for omp_set_lock().
 * @param l Pointer to lock to be set.
 */
inline void autopas_set_lock(autopas_lock_t *l) { omp_set_lock(l); }

/**
 * Wrapper for omp_init_lock().
 * @param l Pointer to lock to be initialized.
 */
inline void autopas_init_lock(autopas_lock_t *l) { omp_init_lock(l); }

/**
 * Wrapper for omp_unset_lock().
 * @param l Pointer to lock to be unset.
 */
inline void autopas_unset_lock(autopas_lock_t *l) { omp_unset_lock(l); }

/**
 * Wrapper for omp_destroy_lock().
 * @param l Pointer to lock to be destroyed.
 */
inline void autopas_destroy_lock(autopas_lock_t *l) { omp_destroy_lock(l); }

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
 * Dummy for omp_lock_t when no OpenMP is available.
 */
typedef int autopas_lock_t;

/**
 * Dummy for omp_set_lock() when no OpenMP is available.
 * @param l Pointer to lock to be set.
 */
inline void autopas_set_lock(autopas_lock_t *l) {
  assert(*l == 0);  // @todo: customize asserts
  *l = 1;
}

/**
 * Dummy for omp_init_lock() when no OpenMP is available.
 * @param l Pointer to lock to be initialized.
 */
inline void autopas_init_lock(autopas_lock_t *l) {
  assert(l != nullptr);  // @todo: customize asserts
}

/**
 * Dummy for omp_unset_lock() when no OpenMP is available.
 * @param l Pointer to lock to be unset.
 */
inline void autopas_unset_lock(autopas_lock_t *l) {
  assert(*l == 1);  // @todo: customize asserts
  *l = 0;
}

/**
 * Dummy for omp_destroy_lock() when no OpenMP is available.
 * @param l Pointer to lock to be destroyed.
 */
inline void autopas_destroy_lock(autopas_lock_t *l) {
  assert(l != nullptr);  // @todo: customize asserts
}

#endif

}  // namespace autopas