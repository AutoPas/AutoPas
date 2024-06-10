/**
 * @file FuzzySetFactory.cpp
 * @author Manuel Lerchner
 * @date 19.04.24
 */

#include "FuzzySetFactory.h"

#include <cmath>
#include <numeric>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::fuzzy_logic {

std::map<std::string, FuzzySetFactory::AvailableFunctions> FuzzySetFactory::availableFunctionMap = {
    {"Triangle", FuzzySetFactory::AvailableFunctions::Triangle},
    {"Trapezoid", FuzzySetFactory::AvailableFunctions::Trapezoid},
    {"Gaussian", FuzzySetFactory::AvailableFunctions::Gaussian},
    {"Sigmoid", FuzzySetFactory::AvailableFunctions::Sigmoid},
    {"SigmoidFinite", FuzzySetFactory::AvailableFunctions::SigmoidFinite}};

std::shared_ptr<FuzzySet> FuzzySetFactory::makeFuzzySet(const std::string &linguisticTerm,
                                                        const std::string &functionName,
                                                        const std::vector<double> &params) {
  // Check if the function name is supported.
  if (availableFunctionMap.find(functionName) == availableFunctionMap.end()) {
    std::string supportedFunctions =
        std::accumulate(availableFunctionMap.begin(), availableFunctionMap.end(), std::string(),
                        [](const std::string &acc, const std::pair<const std::string, AvailableFunctions> &b) {
                          return acc + b.first + ", ";
                        });
    autopas::utils::ExceptionHandler::exception(
        "Cannot create FuzzySet: {}. Unsupported function name: {}. Choose from {}", linguisticTerm, functionName,
        supportedFunctions);
  }

  FuzzySet::BaseMembershipFunction baseMembershipFunction;

  // Create the membership function based on the function name.
  switch (availableFunctionMap.at(functionName)) {
    case AvailableFunctions::Triangle:
      if (params.size() != 3) {
        throwInvalidNumberOfArguments(linguisticTerm, functionName, 3, params.size());
      }
      baseMembershipFunction = std::tuple(functionName, params, triangleFunction(params[0], params[1], params[2]));
      break;
    case AvailableFunctions::Trapezoid:
      if (params.size() != 4) {
        throwInvalidNumberOfArguments(linguisticTerm, functionName, 4, params.size());
      }
      baseMembershipFunction =
          std::tuple(functionName, params, trapezoidFunction(params[0], params[1], params[2], params[3]));
      break;
    case AvailableFunctions::Gaussian:
      if (params.size() != 2) {
        throwInvalidNumberOfArguments(linguisticTerm, functionName, 2, params.size());
      }
      baseMembershipFunction = std::tuple(functionName, params, gaussianFunction(params[0], params[1]));
      break;
    case AvailableFunctions::Sigmoid:
      if (params.size() != 2) {
        throwInvalidNumberOfArguments(linguisticTerm, functionName, 2, params.size());
      }
      baseMembershipFunction = std::tuple(functionName, params, sigmoidFunction(params[0], params[1]));
      break;
    case AvailableFunctions::SigmoidFinite:
      if (params.size() != 3) {
        throwInvalidNumberOfArguments(linguisticTerm, functionName, 3, params.size());
      }
      baseMembershipFunction = std::tuple(functionName, params, sigmoidFiniteFunction(params[0], params[1], params[2]));
      break;
    default:
      autopas::utils::ExceptionHandler::exception("Unknown function name");
  }

  return std::make_shared<FuzzySet>(linguisticTerm, std::move(baseMembershipFunction));
}

std::function<double(double)> FuzzySetFactory::triangleFunction(double min, double mid, double max) {
  auto triangular = [min, mid, max](double value) {
    if (value <= min or value >= max) {
      return 0.0;
    } else if (value <= mid) {
      return (value - min) / (mid - min);
    } else {
      return (max - value) / (max - mid);
    }
  };
  return triangular;
}

std::function<double(double)> FuzzySetFactory::trapezoidFunction(double min, double mid1, double mid2, double max) {
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
  return trapezoidal;
}

std::function<double(double)> FuzzySetFactory::gaussianFunction(double mean, double sigma) {
  auto gaussian = [mean, sigma](double value) { return std::exp(-0.5 * std::pow((value - mean) / sigma, 2)); };
  return gaussian;
}

std::function<double(double)> FuzzySetFactory::sigmoidFunction(double center, double slope) {
  auto sigmoid = [center, slope](double value) { return 1.0 / (1.0 + std::exp(-slope * (value - center))); };
  return sigmoid;
}

std::function<double(double)> FuzzySetFactory::sigmoidFiniteFunction(double lower, double center, double upper) {
  auto s_shape = [](double x, double lower, double center, double upper) {
    if (x <= lower) {
      return 0.0;
    } else if (lower <= x && x <= center) {
      return 0.5 * std::pow((x - lower) / (center - lower), 2);
    } else if (center <= x && x <= upper) {
      return 1.0 - 0.5 * std::pow((x - upper) / (center - upper), 2);
    } else {
      return 1.0;
    }
  };

  auto sigmoidFinite = [lower, center, upper, s_shape](double value) {
    bool is_reversed = lower > upper;

    if (is_reversed) {
      return 1.0 - s_shape(value, upper, center, lower);
    } else {
      return s_shape(value, lower, center, upper);
    }
  };
  return sigmoidFinite;
}

void FuzzySetFactory::throwInvalidNumberOfArguments(const std::string &linguisticTerm, const std::string &functionName,
                                                    size_t expected, size_t actual) {
  autopas::utils::ExceptionHandler::exception(
      "Cannot create FuzzySet: {}. Wrong number of parameters for {}: Expected {}, got {}", linguisticTerm,
      functionName, expected, actual);
}

}  //   namespace autopas::fuzzy_logic