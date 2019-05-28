/**
 * @file ArrayMath.h
 *
 * @date 18 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <numeric>

namespace autopas {

/**
 * Class to handle mathematical operations of floating point std::array's.
 */
class ArrayMath {
 public:
  /**
   * adds two arrays, returns the result.
   * @tparam T floating point type
   * @tparam SIZE size of the arrays
   * @param a first summand
   * @param b second summand
   * @return a + b
   */
  template <class T, std::size_t SIZE>
  static std::array<T, SIZE> add(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    std::array<T, SIZE> result;
    for (std::size_t d = 0; d < SIZE; ++d) {
      result[d] = a[d] + b[d];
    }
    return result;
  }

  /**
   * Subtracts array b from array a and returns the result.
   * @tparam T floating point type
   * @tparam SIZE size of the arrays
   * @param a
   * @param b
   * @return a - b
   */
  template <class T, std::size_t SIZE>
  static std::array<T, SIZE> sub(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    std::array<T, SIZE> result;
    for (std::size_t d = 0; d < SIZE; ++d) {
      result[d] = a[d] - b[d];
    }
    return result;
  }

  /**
   * Multiplies two array's element wise and returns the result.
   * @tparam T floating point type
   * @tparam SIZE size of the arrays
   * @param a
   * @param b
   * @return element-wise multiplication of a and b
   */
  template <class T, std::size_t SIZE>
  static std::array<T, SIZE> mul(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    std::array<T, SIZE> result;
    for (std::size_t d = 0; d < SIZE; ++d) {
      result[d] = a[d] * b[d];
    }
    return result;
  }

  /**
   * Adds a scalar s to each element of array a and returns the result.
   * @tparam T floating point type
   * @tparam SIZE size of the array a
   * @param a the array
   * @param s the scalar to be added to each element of a
   * @return array who's elements are a[i]+s
   */
  template <class T, std::size_t SIZE>
  static std::array<T, SIZE> addScalar(const std::array<T, SIZE> &a, T s) {
    std::array<T, SIZE> result;
    for (std::size_t d = 0; d < SIZE; ++d) {
      result[d] = a[d] + s;
    }
    return result;
  }

  /**
   * Subtracts a scalar s from each element of array a and returns the result.
   * @tparam T floating point type
   * @tparam SIZE size of the array a
   * @param a the array
   * @param s the scalar to be subtracted from each element of a
   * @return array who's elements are a[i]-s
   */
  template <class T, std::size_t SIZE>
  static std::array<T, SIZE> subScalar(const std::array<T, SIZE> &a, T s) {
    std::array<T, SIZE> result;
    for (std::size_t d = 0; d < SIZE; ++d) {
      result[d] = a[d] - s;
    }
    return result;
  }

  /**
   * Multiplies a scalar s to each element of array a and returns the result.
   * @tparam T floating point type
   * @tparam SIZE size of the array a
   * @param a the array
   * @param s the scalar to be multiplied to each element of a
   * @return array who's elements are a[i]*s
   */
  template <class T, std::size_t SIZE>
  static std::array<T, SIZE> mulScalar(const std::array<T, SIZE> &a, T s) {
    std::array<T, SIZE> result;
    for (std::size_t d = 0; d < SIZE; ++d) {
      result[d] = a[d] * s;
    }
    return result;
  }

  /**
   * Generates the dot product of two arrays.
   * Returns the sum of a[i]*b[i] summed over all i, where i is in [0, SIZE)
   * @tparam T floating point type
   * @tparam SIZE size of the arrays
   * @param a first array
   * @param b second array
   * @return dot product of a and b
   */
  template <class T, std::size_t SIZE>
  static T dot(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    auto result = static_cast<T>(0.0);
    for (std::size_t d = 0; d < SIZE; ++d) {
      result += a[d] * b[d];
    }
    return result;
  }

  /**
   * Generates a normalized array (|a| = 1).
   * @tparam T floating point type
   * @tparam SIZE size of the array
   * @param a input array
   * @return normalized array of a
   */
  template <class T, std::size_t SIZE>
  static std::array<T, SIZE> normalize(const std::array<T, SIZE> &a) {
    const T sum = std::accumulate(a.cbegin(), a.cend(), static_cast<T>(0));
    return mulScalar(a, static_cast<T>(1) / sum);
  }

};  // class ArrayMath

}  // namespace autopas
