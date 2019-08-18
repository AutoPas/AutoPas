/**
 * @file Math.h
 * @author Jan Nguyen
 * @date 18.08.19
 */

#pragma once

#include <cmath>

namespace autopas::Math {
/**
 * Probability density function PDF of the standard normal distribution
 * @param x
 * @return PDF(x)
 */
double normalPDF(double x);

/**
 * Cumulative distribution function CDF of the standard normal distribution
 * @param x
 * @return CDF(x)
 */
double normalCDF(double x);

}  // namespace autopas::Math
