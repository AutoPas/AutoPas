/**
 * @file AVXUtils.h
 * @date 17/05/2022
 * @author Q. Behrami
 */
#pragma once

#ifndef __AVX__
#pragma message "AVXUtils.h included, but AVX is not supported by the compiler."
#else
#include "immintrin.h"

namespace autopas::utils::avx {

/**
 * @brief Utility function to calculate the sum of an AVX register horizontally.
 * @param data register of 4 doubles
 * @return sum of its elements
 */
inline double horizontalSum(const __m256d &data) {
  __m256d sum = _mm256_hadd_pd(data, data);
  __m256d permuted = _mm256_permute4x64_pd(sum, _MM_PERM_ACBD);
  sum = _mm256_hadd_pd(sum, permuted);
  sum = _mm256_permute4x64_pd(sum, _MM_PERM_BBDD);
  return _mm_cvtsd_f64(_mm256_extractf128_pd(sum, 0));
}

/**
 * @brief Utility function to conditionally load contiguous 64-bit doubles from memory
 * @param useMask whether to use the mask
 * @param data the data to load from
 * @param mask the mask to use
 * @return the loaded data
 */
inline __m256d load_pd(bool useMask, const double *const __restrict data, const __m256i &mask) {
  if (useMask) {
    return _mm256_maskload_pd(data, mask);
  } else {
    return _mm256_loadu_pd(data);
  }
}

/**
 * @brief Utility function to conditionally load contiguous 64-bit integers from memory
 * @param useMask whether to use the mask
 * @param data the data to load from
 * @param mask the mask to use
 * @return the loaded data
 */
inline __m256i load_epi64(bool useMask, const size_t *const __restrict data, const __m256i &mask) {
  if (useMask) {
    return _mm256_maskload_epi64(reinterpret_cast<const long long int *>(data), mask);
  } else {
    return _mm256_loadu_si256(reinterpret_cast<const __m256i *>(data));
  }
}

/**
 * @brief Utility function to conditionally store contiguous 64-bit integers into a memory location
 * @param useMask whether to use the mask
 * @param data the data to store to
 * @param mask the mask to use
 * @param values the data to be stored
 */
inline void store_epi64(bool useMask, size_t *const __restrict data, const __m256i &mask, const __m256i &values){
  if (useMask) {
    _mm256_maskstore_epi64(reinterpret_cast<long long int *>(data), mask, values);
  } else {
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(data), values);
  }
}

/**
 * @brief Utility function to conditionally load non-contiguous 64-bit doubles from memory
 * @param useMask whether to use the mask
 * @param data the data to load from
 * @param indices the indices to load data from
 * @param mask the mask to use
 * @return the loaded data
 */
inline __m256d gather_pd(bool useMask, const double *const __restrict data, const __m256i &indices,
                                const __m256i &mask) {
  if (useMask) {
    return _mm256_mask_i64gather_pd(_mm256_setzero_pd(), data, indices, _mm256_castsi256_pd(mask), 8);
  } else {
    return _mm256_i64gather_pd(data, indices, 8);
  }
}

/**
 * @brief Utility function to conditionally load non-contiguous 64-bit integers from memory
 * @param useMask whether to use the mask
 * @param data the data to load from
 * @param indices the indices to load data from
 * @param mask the mask to use
 * @return the loaded data
 */
inline __m256i gather_epi64(bool useMask, const size_t *const __restrict data, const __m256i &indices,
                                   const __m256i &mask) {
  if (useMask) {
    return _mm256_mask_i64gather_epi64(_mm256_setzero_si256(), reinterpret_cast<const long long int*>(data), indices, mask, 8);
  } else {
    return _mm256_i64gather_epi64(reinterpret_cast<const long long int*>(data), indices, 8);
  }
}

/**
 * @brief Utility function to conditionally store non-contiguous 64-bit integers in memory
 * @param useMask whether to use the mask
 * @param ptr the location to store the data
 * @param indices the indices to store data at
 * @param mask the mask to use
 * @note May be improved with AVX512 scatter instructions
 */
inline void scatter_pd(bool useMask, double *const __restrict ptr, const __m256i &indices, const __m256i &mask,
                              __m256d &values) {
  alignas(32) double temp[4];
  alignas(32) size_t idx[4];

  _mm256_store_pd(temp, values);
  _mm256_store_si256((__m256i *)idx, indices);

  for (size_t i = 0; i < 4; ++i) {
    if (mask[i]) {  // <- This does not work for all compilers >:(
      ptr[idx[i]] = temp[i];
    } else
      break;
  }
}
}  // namespace autopas::utils::avx
#endif