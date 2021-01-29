/**
 * @file Math.cpp
 * @author Jan Nguyen
 * @date 18.08.19
 */

#include "Math.h"

namespace autopas::utils::Math {

double normalPDF(double x) {
  const double factor = 1. / std::sqrt(2 * M_PI);
  return factor * std::exp(-x * x / 2.);
}

double normalCDF(double x) {
  const double factor = -1. / std::sqrt(2);
  return std::erfc(factor * x) / 2.;
}

double sigmoid(double x) {
  if (x >= 0) {
    return 1. / (1. + std::exp(-x));
  } else {
    // too big negative values may overflow exp(-x)
    // therefore convert 1/(1+exp(-x)) to exp(x)/(1+exp(x))
    double ex = std::exp(x);
    return ex / (1. + ex);
  }
}

Eigen::VectorXd makeVectorXd(const std::vector<double> &elements) {
  return Eigen::Map<const Eigen::VectorXd>(elements.data(), elements.size());
}

Eigen::VectorXi makeVectorXi(const std::vector<int> &elements) {
  return Eigen::Map<const Eigen::VectorXi>(elements.data(), elements.size());
}

}  // namespace autopas::utils::Math
