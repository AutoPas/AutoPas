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

double roundFixed(double d, int fixedPrecision) {
  const auto factor = std::pow(10, fixedPrecision);
  return std::round(d * factor) / factor;
}

double roundFloating(double d, int floatingPrecision) {
  // Special case for 0 to avoid log10(0)
  if (d == 0.0) {
    return d;
  }
  // The shift factor is the precision minus the number of digits d has before the decimal point
  const auto factor = std::pow(10, floatingPrecision - std::ceil(std::log10(std::abs(d))));
  return std::round(d * factor) / factor;
}

Eigen::VectorXd makeVectorXd(const std::vector<double> &elements) {
  return Eigen::Map<const Eigen::VectorXd>(elements.data(), elements.size());
}

Eigen::VectorXi makeVectorXi(const std::vector<int> &elements) {
  return Eigen::Map<const Eigen::VectorXi>(elements.data(), elements.size());
}

}  // namespace autopas::utils::Math
