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

namespace autopas::ArrayMath {

/**
 * Adds two arrays, returns the result.
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a first summand
 * @param b second summand
 * @return a + b
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> add(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
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
constexpr std::array<T, SIZE> sub(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
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
constexpr std::array<T, SIZE> mul(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  std::array<T, SIZE> result{};
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
constexpr std::array<T, SIZE> addScalar(const std::array<T, SIZE> &a, T s) {
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
constexpr std::array<T, SIZE> subScalar(const std::array<T, SIZE> &a, T s) {
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
constexpr std::array<T, SIZE> mulScalar(const std::array<T, SIZE> &a, T s) {
  std::array<T, SIZE> result{};
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
constexpr T dot(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  auto result = static_cast<T>(0.0);
  for (std::size_t d = 0; d < SIZE; ++d) {
    result += a[d] * b[d];
  }
  return result;
}

/**
 * Creates a new array by performing an element-wise static_cast<>.
 * @tparam output_t Output type.
 * @tparam input_t Input type.
 * @tparam SIZE Size of the array.
 * @param a Input array.
 * @return Array of type std::array<output_t, SIZE>.
 */
template <class output_t, class input_t, std::size_t SIZE>
constexpr std::array<output_t, SIZE> static_cast_array(const std::array<input_t, SIZE> &a) {
  std::array<output_t, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = static_cast<output_t>(a[d]);
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
constexpr std::array<T, SIZE> normalize(const std::array<T, SIZE> &a) {
  const T length = std::sqrt(dot(a, a));
  return mulScalar(a, static_cast<T>(1) / length);
}

/**
 * Generates a string representation of a container which fulfills the Container requirement (provide cbegin and cend).
 * @tparam T Type of Container.
 * @param a Container.
 * @param delimiter String delimiter.
 * @return String representation of a.
 */
template <class T>
std::string to_string(const T &a, const std::string &delimiter = ", ") {
  auto it = std::cbegin(a);
  const auto end = std::cend(a);
  if (it == end) {
    return "";
  }
  std::ostringstream strStream;
  strStream << *it;
  for (++it; it != end; ++it) {
    strStream << delimiter << *it;
  }

  return strStream.str();
}

}  // namespace autopas::ArrayMath
