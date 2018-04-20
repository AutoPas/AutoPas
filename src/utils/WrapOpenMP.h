/**
 * @file WrapOpenMP.h
 * @author F. Gratl
 * @date 4/20/18
 */

#pragma once

/**
 * Provide non-OpenMP versions of the most common OpenMP function calls,
 * so that they don't have to be wrapped in #ifdef-s every time.
 *
 * Proper wrapper and renaming necessary, because of -fopenmp-simd handling of
 * gcc.
 *
 * Extend when necessary.
 */

#include <cassert>
#if defined(_OPENMP)
#include <omp.h>

inline int autopas_get_thread_num() { return omp_get_thread_num(); }
inline int autopas_get_num_threads() { return omp_get_num_threads(); }
inline int autopas_get_max_threads() { return omp_get_max_threads(); }
typedef omp_lock_t autopas_lock_t;
inline void autopas_set_lock(autopas_lock_t* l) { omp_set_lock(l); }
inline void autopas_init_lock(autopas_lock_t* l) { omp_init_lock(l); }
inline void autopas_unset_lock(autopas_lock_t* l) { omp_unset_lock(l); }
inline void autopas_destroy_lock(autopas_lock_t* l) { omp_destroy_lock(l); }

#else
inline int autopas_get_thread_num() { return 0; }
inline int autopas_get_num_threads() { return 1; }
inline int autopas_get_max_threads() { return 1; }
typedef int autopas_lock_t;
inline void autopas_set_lock(autopas_lock_t *l) {
  assert(*l == 0); // TODO: customize asserts
  *l = 1;
}
inline void autopas_init_lock(autopas_lock_t *l) {
  assert(l != nullptr); // TODO: customize asserts
}
inline void autopas_unset_lock(autopas_lock_t *l) {
  assert(*l == 1); // TODO: customize asserts
  *l = 0;
}
inline void autopas_destroy_lock(autopas_lock_t *l) {
  assert(l != nullptr); // TODO: customize asserts
}
#endif