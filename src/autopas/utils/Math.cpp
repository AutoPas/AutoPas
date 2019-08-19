/**
 * @file Math.cpp
 * @author Jan Nguyen
 * @date 18.08.19
 */

#include "Math.h"

namespace autopas::Math {
/**
 * Probability density function PDF of the standard normal distribution
 * @param x
 * @return PDF(x)
 */
double normalPDF(double x) {
  static double factor = 1. / std::sqrt(2 * M_PI);
  return factor * std::exp(-x * x / 2.);
}

/**
 * Cumulative distribution function CDF of the standard normal distribution
 * @param x
 * @return CDF(x)
 */
double normalCDF(double x) {
  static double factor = -1. / std::sqrt(2);
  return std::erfc(factor * x) / 2.;
}

}  // namespace autopas::Math
