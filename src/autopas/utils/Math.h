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
