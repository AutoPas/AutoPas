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
 * True iff for all d the following holds: a[d] < b[d]
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a first parameter
 * @param b second parameter
 * @return a < b
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr bool less(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  bool result = true;
  for (std::size_t d = 0; d < SIZE; ++d) {
    result = result and (a[d] < b[d]);
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
  T result = 0.;
  for (std::size_t i = 0; i < SIZE; i++) {
    result += a[i] * b[i];
  }
  return result;
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

/**
 * For each element in a, computes the smallest integer value not less than the element.
 * @tparam T floating point type
 * @tparam SIZE size of the array
 * @param a input array
 * @return rounded up values of a
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<T, SIZE> ceil(const std::array<T, SIZE> &a) {
  std::array<T, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = std::ceil(a[d]);
  }
  return result;
}

/**
 * Floors all array elements and converts them to integers.
 * @tparam T floating point type
 * @tparam SIZE size of the array
 * @param a input array
 * @return New array with floored elements of new type int.
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<int, SIZE> floorToInt(const std::array<T, SIZE> &a) {
  std::array<int, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = static_cast<int>(std::floor(a[d]));
  }
  return result;
}

/**
 * Ceils all array elements and converts them to integers.
 * @tparam T floating point type
 * @tparam SIZE size of the array
 * @param a input array
 * @return New array with ceiled elements of new type int.
 */
template <class T, std::size_t SIZE>
[[nodiscard]] constexpr std::array<int, SIZE> ceilToInt(const std::array<T, SIZE> &a) {
  std::array<int, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = static_cast<int>(std::ceil(a[d]));
  }
  return result;
}

// namespace for templated operators
inline namespace literals {

/**
 * Adds two arrays, returns the result.
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a first summand
 * @param b second summand
 * @return a + b
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> operator+(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  return add(a, b);
}

/**
 * Assignment operator to add two arrays
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a first summand
 * @param b second summand
 * @return a + b
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> &operator+=(std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  for (std::size_t d = 0; d < SIZE; ++d) {
    a[d] += b[d];
  }
  return a;
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
constexpr std::array<T, SIZE> operator-(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  return sub(a, b);
}

/**
 * Assignment operator to subtract two arrays
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a
 * @param b
 * @return a - b
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> &operator-=(std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  for (std::size_t d = 0; d < SIZE; ++d) {
    a[d] -= b[d];
  }
  return a;
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
constexpr std::array<T, SIZE> operator*(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  return mul(a, b);
}

/**
 * Assignment operator to multply two arrays
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a
 * @param b
 * @return element-wise multiplication of a and b
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> &operator*=(std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  for (std::size_t d = 0; d < SIZE; ++d) {
    a[d] *= b[d];
  }
  return a;
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
constexpr std::array<T, SIZE> operator/(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  return div(a, b);
}

/**
 * Assignment operator to divide two arrays
 * @tparam T floating point type
 * @tparam SIZE size of the arrays
 * @param a dividend.
 * @param b divisor.
 * @return element-wise quotient of a and b, i.e., `result[i] = a[i]/b[i]`
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> &operator/=(std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
  for (std::size_t d = 0; d < SIZE; ++d) {
    a[d] /= b[d];
  }
  return a;
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
constexpr std::array<T, SIZE> operator+(const std::array<T, SIZE> &a, T s) {
  return addScalar(a, s);
}

/**
 * Assignment operator to add a scalar s to each element of array
 * @tparam T floating point type
 * @tparam SIZE size of the array a
 * @param a the array
 * @param s the scalar to be added to each element of a
 * @return array who's elements are a[i]+s
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> &operator+=(std::array<T, SIZE> &a, T s) {
  for (std::size_t d = 0; d < SIZE; ++d) {
    a[d] += s;
  }
  return a;
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
constexpr std::array<T, SIZE> operator-(const std::array<T, SIZE> &a, T s) {
  return subScalar(a, s);
}

/**
 * Assignment operator to subtract a scalar s to each element of array
 * @tparam T floating point type
 * @tparam SIZE size of the array a
 * @param a the array
 * @param s the scalar to be subtracted from each element of a
 * @return array who's elements are a[i]-s
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> &operator-=(std::array<T, SIZE> &a, T s) {
  for (std::size_t d = 0; d < SIZE; ++d) {
    a[d] -= s;
  }
  return a;
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
constexpr std::array<T, SIZE> operator*(const std::array<T, SIZE> &a, T s) {
  return mulScalar(a, s);
}

/**
 * Assignment operator to multiply a scalar s to each element of array
 * @tparam T floating point type
 * @tparam SIZE size of the array a
 * @param a the array
 * @param s the scalar to be multiplied to each element of a
 * @return array who's elements are a[i]*s
 */
template <class T, std::size_t SIZE>
constexpr std::array<T, SIZE> &operator*=(std::array<T, SIZE> &a, T s) {
  for (std::size_t d = 0; d < SIZE; ++d) {
    a[d] *= s;
  }
  return a;
}

}  // namespace literals

}  // namespace autopas::utils::ArrayMath
