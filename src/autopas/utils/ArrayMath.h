/**
 * @file ArrayMath.h
 *
 * @date 18 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <sstream>

namespace autopas::utils::ArrayMath {

/**
 * Adds two arrays, returns the result.
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a first summand
 * @param b second summand
 * @return a + b
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<T, SIZE> add(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
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
[[nodiscard]] constexpr std::array<T, SIZE> sub(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = a[d] - b[d];
  }
  return result;
}

/**
 * Takes elementwise minimum, returns the result.
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a first parameter
 * @param b second parameter
 * @return min(a, b)
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<T, SIZE> min(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = std::min<T>(a[d], b[d]);
  }
  return result;
}

/**
 * Takes elementwise maximum and returns the result.
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a
 * @param b
 * @return max(a, b)
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<T, SIZE> max(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = std::max<T>(a[d], b[d]);
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
[[nodiscard]] constexpr std::array<T, SIZE> mul(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = a[d] * b[d];
  }
  return result;
}

/**
 * Divides two array's element-wise and returns the result.
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a dividend.
 * @param b divisor.
 * @return element-wise quotient of a and b, i.e., `result[i] = a[i]/b[i]`
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<T, SIZE> div(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = a[d] / b[d];
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
[[nodiscard]] constexpr std::array<T, SIZE> addScalar(const std::array<T, SIZE> &a, T s) {
  std::array<T, SIZE> result{};
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
[[nodiscard]] constexpr std::array<T, SIZE> subScalar(const std::array<T, SIZE> &a, T s) {
  std::array<T, SIZE> result{};
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
[[nodiscard]] constexpr std::array<T, SIZE> mulScalar(const std::array<T, SIZE> &a, T s) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = a[d] * s;
  }
  return result;
}

// unnamed namespace for private helper functions
namespace {
/**
 * Helper function to provide a templated dot product that is basically the same as writing it out by hand.
 * @tparam T
 * @tparam I
 * @param a
 * @param b
 * @return
 */
template <typename T, size_t... I>
double dotAux(T a, T b, std::integer_sequence<size_t, I...>) {
  return ((std::get<I>(a) * std::get<I>(b)) + ...);
}
}  // namespace

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
[[nodiscard]] constexpr T dot(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  return dotAux(a, b, std::make_index_sequence<SIZE>{});
}

/**
 * Generates the cross product of two arrays of 3 floats.
 * @tparam T floating point type
 * @param a 3D vector (denoted by array of 3 floats)
 * @param b 3D vector (denoted by array of 3 floats)
 * @return cross product a x b
 */
template <class T>
[[nodiscard]] constexpr std::array<T, 3> cross(const std::array<T, 3> &a, const std::array<T, 3> &b) {
  return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}

/**
 * Computes the L2Norm / Euclidean norm.
 * @tparam T floating point type
 * @tparam SIZE size of the array
 * @param a input array
 * @return L2Norm
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr T L2Norm(const std::array<T, SIZE> &a) {
  return std::sqrt(dot(a, a));
}

/**
 * Computes the product of all elements in a.
 * @tparam T input array type, which fulfills the C++ container requirement
 * @param a input array
 * @return product
 */
template <class T>
[[nodiscard]] constexpr typename T::value_type prod(const T &a) {
  return std::accumulate(a.cbegin(), a.cend(), static_cast<typename T::value_type>(1), std::multiplies<>());
}

/**
 * Computes the absolute value for all elements in a.
 * @tparam T floating point type
 * @tparam SIZE size of the array
 * @param a input array
 * @return absolute values of a
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<T, SIZE> abs(const std::array<T, SIZE> &a) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = std::abs(a[d]);
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
[[nodiscard]] constexpr std::array<T, SIZE> normalize(const std::array<T, SIZE> &a) {
  return mulScalar(a, static_cast<T>(1) / L2Norm(a));
}

}  // namespace autopas::utils::ArrayMath
