/**
 * @file ConstexprMath.h
 *
 * @date 15.05.2023
 * @author D. Martin
 */

#pragma once

#include <limits>

namespace autopas::utils::ConstexprMath {

/**
 * Calculates the square root of floating point values x based on Newton-Raphson methon
 * @tparam T floating point type
 * @param x input value
 * @param epsilon epsilon value used for floating point accuracy comparison
 * @return sqrt(x)
 */
template <typename T>
constexpr typename std::enable_if_t<std::is_floating_point_v<T>, T> sqrt(T x, T epsilon) {
  if (x >= 0 and x < std::numeric_limits<T>::infinity()) {
    T xn = x;
    T prev = 0;
    // the while loop checks if the absolut value of (xn - prev) is greater as the specified epsilon
    while ((xn - prev < 0 ? prev - xn : xn - prev) > epsilon) {
      prev = xn;
      T xn1 = 0.5 * (xn + x / xn);
      xn = xn1;
    }
    return xn;
  } else {
    return std::numeric_limits<T>::quiet_NaN();
  }
}

/**
 * Calculates the square root of integral values x based on Newton-Raphson methon
 * @tparam T integral type
 * @param x input value
 * @return sqrt(x)
 */
template <typename T>
constexpr typename std::enable_if_t<std::is_integral_v<T>, T> sqrt(T x) {
  // see https://en.wikipedia.org/wiki/Integer_square_root#Example_implementation_in_C
  if (x >= 0) {
    if (x <= 1) {
      return x;
    }
    T x0 = x / 2;
    T x1 = (x0 + x / x0) / 2;
    while (x1 < x0) {
      x0 = x1;
      x1 = (x0 + x / x0) / 2;
    }
    return x0;
  } else {
    throw std::invalid_argument("Negative number passed.");
  }
}

}  // namespace autopas::utils::ConstexprMath