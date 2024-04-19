/**
 * @file ComposedMembershipFunction.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "MembershipFunction.h"

#include <cmath>

namespace autopas::fuzzy_logic {

std::shared_ptr<FuzzySet> makeTriangle(const std::string &linguisticTerm, double min, double mid, double max) {
  auto triangular = std::make_shared<FuzzySet::BaseMembershipFunction>([min, mid, max](double value) {
    if (value <= min or value >= max) {
      return 0.0;
    } else if (value <= mid) {
      return (value - min) / (mid - min);
    } else {
      return (max - value) / (max - mid);
    }
  });
  return std::make_shared<FuzzySet>(linguisticTerm, triangular);
}

std::shared_ptr<FuzzySet> makeTrapezoid(const std::string &linguisticTerm, double min, double mid1, double mid2,
                                        double max) {
  auto trapezoidal = std::make_shared<FuzzySet::BaseMembershipFunction>([min, mid1, mid2, max](double value) {
    if (value <= min or value >= max) {
      return 0.0;
    } else if (value <= mid1) {
      return (value - min) / (mid1 - min);
    } else if (value <= mid2) {
      return 1.0;
    } else {
      return (max - value) / (max - mid2);
    }
  });
  return std::make_shared<FuzzySet>(linguisticTerm, trapezoidal);
}

std::shared_ptr<FuzzySet> makeGaussian(const std::string &linguisticTerm, double mean, double sigma) {
  auto gaussian = std::make_shared<FuzzySet::BaseMembershipFunction>(
      [mean, sigma](double value) { return std::exp(-0.5 * std::pow((value - mean) / sigma, 2)); });
  return std::make_shared<FuzzySet>(linguisticTerm, gaussian);
}

}  //   namespace autopas::fuzzy_logic