/**
 * @file FuzzySetFactory.cpp
 * @author Manuel Lerchner
 * @date 19.04.24
 */

#include <cmath>

#include "FuzzySetFactory.h"

namespace autopas::fuzzy_logic {

std::shared_ptr<FuzzySet> makeTriangle(const std::string &linguisticTerm, double min, double mid, double max) {
  auto triangular = [min, mid, max](double value) {
    if (value <= min or value >= max) {
      return 0.0;
    } else if (value <= mid) {
      return (value - min) / (mid - min);
    } else {
      return (max - value) / (max - mid);
    }
  };
  return std::make_shared<FuzzySet>(linguisticTerm, triangular);
}

std::shared_ptr<FuzzySet> makeTrapezoid(const std::string &linguisticTerm, double min, double mid1, double mid2,
                                        double max) {
  auto trapezoidal = [min, mid1, mid2, max](double value) {
    if (value <= min or value >= max) {
      return 0.0;
    } else if (value <= mid1) {
      return (value - min) / (mid1 - min);
    } else if (value <= mid2) {
      return 1.0;
    } else {
      return (max - value) / (max - mid2);
    }
  };
  return std::make_shared<FuzzySet>(linguisticTerm, trapezoidal);
}

std::shared_ptr<FuzzySet> makeGaussian(const std::string &linguisticTerm, double mean, double sigma) {
  auto gaussian = [mean, sigma](double value) { return std::exp(-0.5 * std::pow((value - mean) / sigma, 2)); };
  return std::make_shared<FuzzySet>(linguisticTerm, gaussian);
}

}  //   namespace autopas::fuzzy_logic