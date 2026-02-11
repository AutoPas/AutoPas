/**
 * @file FloatComparison.h
 * @author J. Schuhmacher
 * @date 11.02.26
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <type_traits>

namespace autopas::utils {

/**
 * This default relative EPSILON used for the {@link almostEqualRelative} function.
 */
constexpr double EPSILON_ALMOST_EQUAL = 1e-10;

/**
 * The default maximal allowed ULP distance utilized for FloatingPoint comparisons using the
 * {@link almostEqualUlps} function. This is also the value utilized in GoogleTest for their DoubleEq() and FloatEq()
 * Matchers.
 *
 * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */
constexpr unsigned int MAX_ULP_DISTANCE = 4;

namespace internal {
/**
 * Helper to find out the bit size of a given type
 * @tparam T any type T
 */
template <typename T>
constexpr size_t bit_size_v = sizeof(T) * 8;

/**
 * Helper struct to get the correspondingly sized integer for a given floating type.
 * As 128-bit type floating points are not yet available in C++20, we ensure it will fail to be instantiated - so that
 * the necessary code changes in the conditional can happen.
 * @tparam FloatType floating point type
 */
template <std::floating_point FloatType>
struct int_t_impl {
  static_assert(bit_size_v<FloatType> == 32 || bit_size_v<FloatType> == 64,
                "int_t only supports 32-bit and 64-bit floating-point types.");
  using type = std::conditional_t<bit_size_v<FloatType> == 32, std::int32_t, std::int64_t>;
};

}  // namespace internal

/**
 * Returns the correspondgly sized integer type for a given float. E.g. int_t of float would be int32_t.
 */
template <std::floating_point FloatType>
using int_t = typename internal::int_t_impl<FloatType>::type;

/**
 * Function for comparing closeness of two floating point numbers using ULP (Units in the Last Place) method.
 *
 * @tparam FloatType must be either double or float (ensured by static assertion)
 * @param lhs The left hand side floating point number to compare.
 * @param rhs The right hand side floating point number to compare.
 * @param ulpDistance The maximum acceptable ULP distance between the two floating points
 *      for which they would be considered near each other. This is optional and by default, it will be {@link
 * MAX_ULP_DISTANCE}.
 *
 * @return true if the ULP distance between lhs and rhs is less than or equal to the provided ulpDistance value,
 * otherwise, false. Returns true if both numbers are exactly the same. Returns false if the signs do not match.
 * @note The ULP distance between 3.0 and std::nextafter(3.0, INFINITY) would be 1,
 *      the ULP distance of 3.0 and std::nextafter(std::nextafter(3.0, INFINITY), INFINITY) would be 2, etc.
 * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */
template <std::floating_point FloatType>
bool almostEqualUlps(FloatType lhs, FloatType rhs, unsigned int ulpDistance = MAX_ULP_DISTANCE) {
  // In case the floats are equal in their bitwise representation, return true
  if (lhs == rhs) {
    return true;
  }

  // In case the signs mismatch, return false
  if (lhs < static_cast<FloatType>(0.0) && rhs > static_cast<FloatType>(0.0) ||
      lhs > static_cast<FloatType>(0.0) && rhs < static_cast<FloatType>(0.0)) {
    return false;
  }

  // Compute ULP distance by interpreting the floating point as an equivalent-sized integer
  return reinterpret_cast<int_t<FloatType> &>(rhs) - reinterpret_cast<int_t<FloatType> &>(lhs) <= ulpDistance;
}

/**
 * Function to check if two floating point numbers are relatively equal to each other within a given error range or
 * tolerance.
 *
 * @tparam FloatType must be either double or float (ensured by static assertion)
 * @param lhs The first floating-point number to be compared.
 * @param rhs The second floating-point number to be compared.
 * @param epsilon The tolerance for comparison. Two numbers that are less than epsilon apart are considered equal.
 *                The default value is {@link EPSILON_ALMOST_EQUAL}.
 *
 * @return boolean value - Returns `true` if the absolute difference between `lhs` and `rhs` is less than or equal to
 *                         the relative error factored by the larger of the magnitude of `lhs` and `rhs`. Otherwise,
 * `false`.
 * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */
template <std::floating_point FloatType>
bool almostEqualRelative(FloatType lhs, FloatType rhs, FloatType epsilon = EPSILON_ALMOST_EQUAL) {
  const FloatType diff = std::abs(rhs - lhs);
  const FloatType largerValue = std::max(std::abs(rhs), std::abs(lhs));
  return diff <= largerValue * epsilon;
}

}  // namespace autopas::utils