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
 * @tparam T
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Sum or valOverflow or valUnderflow.
 */
template <class T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
T safeAdd(const T &a, const T &b, const T &valUnderflow = std::numeric_limits<T>::min(),
          const T &valOverflow = std::numeric_limits<T>::max()) {
  T result;
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
 * @tparam T
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Returns value in case of overflow.
 * @return Sum or valUnderflow or valOverflow.
 */
template <class T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
T safeAdd(const T &a, const T &b, const T &valUnderflow = -std::numeric_limits<T>::max(),
          const T &valOverflow = std::numeric_limits<T>::max()) {
  T result = a + b;
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
 * @tparam T
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Difference or valOverflow or valUnderflow.
 */
template <class T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
T safeSub(const T &a, const T &b, const T &valUnderflow = std::numeric_limits<T>::min(),
          const T &valOverflow = std::numeric_limits<T>::max()) {
  T result;
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
 * @tparam T
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Returns value in case of overflow.
 * @return Sum or valUnderflow or valOverflow.
 */
template <class T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
T safeSub(const T &a, const T &b, const T &valUnderflow = -std::numeric_limits<T>::max(),
          const T &valOverflow = std::numeric_limits<T>::max()) {
  T result = a - b;
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
 * @tparam T
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Return value in case of overflow.
 * @return Product or valUnderflow or valOverflow.
 */
template <class T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
T safeMul(const T &a, const T &b, const T &valUnderflow = std::numeric_limits<T>::min(),
          const T &valOverflow = std::numeric_limits<T>::max()) {
  T result;
  bool overflow = __builtin_mul_overflow(a, b, &result);
  if (overflow) {
    // if both args are negative this is an underflow.
    if (a < 0 xor b < 0) {
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
 * @tparam T
 * @param a
 * @param b
 * @param valUnderflow Return value in case of underflow.
 * @param valOverflow Returns value in case of overflow.
 * @return Product or valUnderflow or valOverflow.
 */
template <class T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
T safeMul(const T &a, const T &b, const T &valUnderflow = -std::numeric_limits<T>::max(),
          const T &valOverflow = std::numeric_limits<T>::max()) {
  T result = a * b;
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

}  // namespace autopas::utils::Math
