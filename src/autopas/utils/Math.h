/**
 * @file Math.h
 * @author Jan Nguyen
 * @date 18.08.19
 */

#pragma once

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace autopas::utils::Math {
/**
 * Factor of PDF of standard normal distribution.
 */
const double normalScale = 1. / std::sqrt(2 * M_PI);

/**
 * Addition function for integer types that is safe against over and underflow.
 * If over or underflow is detected, the function returns the specified values.
 * @tparam Int_t
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Sum or valOverflow or valUnderflow.
 */
template <class Int_t, std::enable_if_t<std::is_integral_v<Int_t>, bool> = true>
Int_t safeAdd(const Int_t &a, const Int_t &b, const Int_t &valUnderflow = std::numeric_limits<Int_t>::min(),
              const Int_t &valOverflow = std::numeric_limits<Int_t>::max()) {
  Int_t result;
  bool overflow = __builtin_add_overflow(a, b, &result);
  if (overflow) {
    // if both args are negative this is an underflow.
    if (a < 0 and b < 0) {
      result = valUnderflow;
    } else {
      result = valOverflow;
    }
  }
  return result;
}

/**
 * Addition function for floating point types that is safe against over and underflow.
 * If over or underflow is detected, the function returns the specified values.
 *
 * @note Underflow here refers to a value more negative than representable and not the underflow gap around zero.
 *
 * @tparam Float_t
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Sum or valUnderflow or valOverflow.
 */
template <class Float_t, std::enable_if_t<std::is_floating_point_v<Float_t>, bool> = true>
Float_t safeAdd(const Float_t &a, const Float_t &b, const Float_t &valUnderflow = -std::numeric_limits<Float_t>::max(),
                const Float_t &valOverflow = std::numeric_limits<Float_t>::max()) {
  Float_t result = a + b;
  if (std::isinf(result)) {
    if (result > 0) {
      result = valOverflow;
    } else {
      result = valUnderflow;
    }
  }
  return result;
}

/**
 * Subtraction function for integer types that is safe against over and underflow.
 * If over or underflow is detected, the function returns the specified values.
 * @tparam Int_t
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Difference or valOverflow or valUnderflow.
 */
template <class Int_t, std::enable_if_t<std::is_integral_v<Int_t>, bool> = true>
Int_t safeSub(const Int_t &a, const Int_t &b, const Int_t &valUnderflow = std::numeric_limits<Int_t>::min(),
              const Int_t &valOverflow = std::numeric_limits<Int_t>::max()) {
  Int_t result;
  bool overflow = __builtin_sub_overflow(a, b, &result);
  if (overflow) {
    // underflow is only possible if we subtract a positive number. Everything else is an overflow.
    if (b > 0) {
      result = valUnderflow;
    } else {
      result = valOverflow;
    }
  }
  return result;
}

/**
 * Addition function for floating point types that is safe against over and underflow.
 * If over or underflow is detected, the function returns the specified values.
 *
 * @note Underflow here refers to a value more negative than representable and not the underflow gap around zero.
 *
 * @tparam Float_t
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Sum or valUnderflow or valOverflow.
 */
template <class Float_t, std::enable_if_t<std::is_floating_point_v<Float_t>, bool> = true>
Float_t safeSub(const Float_t &a, const Float_t &b, const Float_t &valUnderflow = -std::numeric_limits<Float_t>::max(),
                const Float_t &valOverflow = std::numeric_limits<Float_t>::max()) {
  Float_t result = a - b;
  if (std::isinf(result)) {
    if (result > 0) {
      result = valOverflow;
    } else {
      result = valUnderflow;
    }
  }
  return result;
}

/**
 * Multiplication function for integer types that is safe against over and underflow.
 * If over or underflow is detected, the function returns the specified values.
 * @tparam Int_t
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Product or valUnderflow or valOverflow.
 */
template <class Int_t, std::enable_if_t<std::is_integral_v<Int_t>, bool> = true>
Int_t safeMul(const Int_t &a, const Int_t &b, const Int_t &valUnderflow = std::numeric_limits<Int_t>::min(),
              const Int_t &valOverflow = std::numeric_limits<Int_t>::max()) {
  Int_t result;
  bool overflow = __builtin_mul_overflow(a, b, &result);
  if (overflow) {
    // if exactly one arg is negative this is an underflow.
    if ((a < 0) xor (b < 0)) {
      result = valUnderflow;
    } else {
      result = valOverflow;
    }
  }
  return result;
}

/**
 * Multiplication function for floating point types that is safe against over and underflow.
 * If over or underflow is detected, the function returns the specified values.
 *
 * @note Underflow here refers to a value more negative than representable and not the underflow gap around zero.
 *
 * @tparam Float_t
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Product or valUnderflow or valOverflow.
 */
template <class Float_t, std::enable_if_t<std::is_floating_point_v<Float_t>, bool> = true>
Float_t safeMul(const Float_t &a, const Float_t &b, const Float_t &valUnderflow = -std::numeric_limits<Float_t>::max(),
                const Float_t &valOverflow = std::numeric_limits<Float_t>::max()) {
  Float_t result = a * b;
  if (std::isinf(result)) {
    if (result > 0) {
      result = valOverflow;
    } else {
      result = valUnderflow;
    }
  }
  return result;
}

/**
 * No-overhead power function with exponent known at compile time.
 * @tparam exponent
 * @tparam T
 * @param base
 * @return
 */
template <size_t exponent, class T>
T pow(const T &base) {
  if (exponent == 0) {
    return 1;
  }

  T res = base;
  // the compiler should unroll this loop
  for (size_t i = 0; i < exponent - 1; ++i) {
    res *= base;
  }
  return res;
}

/**
 * Probability density function PDF of the standard normal distribution.
 * @param x
 * @return PDF(x)
 */
double normalPDF(double x);

/**
 * Cumulative distribution function CDF of the standard normal distribution.
 * @param x
 * @return CDF(x)
 */
double normalCDF(double x);

/**
 * Sigmoid logistic function.
 * @param x
 * @return S(x)
 */
double sigmoid(double x);

/**
 * This default relative EPSILON used for the {@link isNearRel} function.
 */
constexpr double EPSILON_RELATIVE_EQUALITY = 1e-9;

/**
 * The default maximal allowed ULP (Units in the Last Place) distance utilized for FloatingPoint comparisons.
 *
 * This is also the value utilized in GoogleTest for their DoubleEq() and FloatEq() Matchers.
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
  /** access the int data type */
  using type = std::conditional_t<bit_size_v<FloatType> == 32, std::int32_t, std::int64_t>;
};
} // namespace internal

/**
 * Returns the correspondingly sized integer type for a given float. E.g. int_t of float would be int32_t.
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
 *      for which they would be considered near each other.
 *
 * @return true if the ULP distance between lhs and rhs is less than or equal to the provided ulpDistance value,
 * otherwise, false. Returns true if both numbers are exactly the same. Returns false if the signs do not match.
 * @note The ULP distance between 3.0 and std::nextafter(3.0, INFINITY) would be 1,
 *      the ULP distance of 3.0 and std::nextafter(std::nextafter(3.0, INFINITY), INFINITY) would be 2, etc.
 * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */
template <std::floating_point FloatType>
bool isInUlp(FloatType lhs, FloatType rhs, unsigned int ulpDistance = MAX_ULP_DISTANCE) {
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
 * Determines if two doubles are near each other. This function should be preferred to comparing with ==.
 * @param a
 * @param b
 * @param maxRelativeDifference inclusive, relative to max(|a|, |b|).
 * @return
 */
template <std::floating_point FloatType>
bool isNearRel(FloatType a, FloatType b, double maxRelativeDifference = EPSILON_RELATIVE_EQUALITY) {
  const auto greaterNumber = std::max(std::abs(a), std::abs(b));
  const auto absoluteDifference = maxRelativeDifference * greaterNumber;
  const auto diff = std::abs(a - b);
  return diff <= absoluteDifference;
}

/**
 * Determines if two doubles are near each other. This function should be preferred to comparing with ==.
 * @param a
 * @param b
 * @param maxAbsoluteDifference inclusive
 * @return
 */
template <std::floating_point FloatType>
bool isNearAbs(FloatType a, FloatType b, double maxAbsoluteDifference) {
  return std::abs(a - b) <= maxAbsoluteDifference;
}

/**
 * Round a floating point number to a given number of decimal digits.
 * @param d Number to round.
 * @param fixedPrecision Number of decimal digits. Negative values lead to rounding of digits left of the decimal.
 * @return d rounded to the given number of digits.
 */
double roundFixed(double d, int fixedPrecision);

/**
 * Round a floating point number to a given floating precision.
 * @param d Number to round.
 * @param floatingPrecision Number of significant digits. Values <0 will return in 0.
 * @return d rounded to the given number of digits.
 */
double roundFloating(double d, int floatingPrecision);

/**
 * Create a vector of doubles from given elements
 * @param elements
 * @return
 */
Eigen::VectorXd makeVectorXd(const std::vector<double> &elements);

/**
 * Create a vector of ints from given elements
 * @param elements
 * @return
 */
Eigen::VectorXi makeVectorXi(const std::vector<int> &elements);
} // namespace autopas::utils::Math
