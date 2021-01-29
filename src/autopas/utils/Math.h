/**
 * @file Math.h
 * @author Jan Nguyen
 * @date 18.08.19
 */

#pragma once

#include <Eigen/Core>
#include <cmath>

namespace autopas::utils::Math {

/**
 * PI
 * @todo c++20: replace with std::numbers::pi, see https://en.cppreference.com/w/cpp/numeric/constants
 */
double getPI() { return std::atan(1.) / 4; };

/**
 * Factor of PDF of standard normal distribution.
 */
const double normalScale = 1. / std::sqrt(2 * getPI());

/**
 * No-overhead power function with exponent known at compile time.
 * @tparam exponent
 * @tparam T
 * @param base
 * @return
 */
template <size_t exponent, class T>
T pow(T base) {
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
