/**
 * @file FuzzySetFactory.h
 * @author Manuel Lerchner
 * @date 19.04.24
 */

#pragma once

#include <functional>
#include <set>

#include "FuzzySet.h"

namespace autopas::fuzzy_logic {

/**
 * A factory class that creates FuzzySets based on the given activation function.
 */
class FuzzySetFactory {
 public:
  /**
   * Deleted constructor to prevent instantiation of this class.
   */
  FuzzySetFactory() = delete;

  /**
   * Enum class that represents the available functions for the FuzzySetFactory.
   *
   * If a new function is added, it must be added here and in the `availableFunctionMap` below.
   */
  enum class AvailableFunctions { Triangle, Trapezoid, Gaussian, Sigmoid, SigmoidFinite };

  /**
   * A map of the form {function_name: function_enum}. Used to map function names to function enums.
   */
  static std::map<std::string, AvailableFunctions> availableFunctionMap;

  /**
   * Constructs a FuzzySet with the given linguistic term, based on the given activation function.
   * @param linguisticTerm The linguistic term of this FuzzySet.
   * @param functionName The name of the function to create.
   * @param params The parameters of the function.
   * @return A shared pointer to the created FuzzySet.
   */
  static std::shared_ptr<FuzzySet> makeFuzzySet(const std::string &linguisticTerm, const std::string &functionName,
                                                std::vector<double> &params);

 private:
  /**
   * Returns a function object that represents a triangular membership function.
   * @param min The left border of the triangle. The membership function will be 0 for values smaller than min.
   * @param mid The peak of the triangle. The membership function will be 1 for values equal to mid.
   * @param max The right border of the triangle. The membership function will be 0 for values greater than max.
   */
  static std::function<double(double)> triangleFunction(double min, double mid, double max);

  /**
   * Returns a function object that represents a trapezoidal membership function.
   * @param min The left border of the trapezoid. The membership function will be 0 for values smaller than min.
   * @param mid1 The left peak of the trapezoid. The membership function will be 1 for values equal to mid1.
   * @param mid2 The right peak of the trapezoid. The membership function will be 1 for values equal to mid2.
   * @param max The right border of the trapezoid. The membership function will be 0 for values greater than max.
   */
  static std::function<double(double)> trapezoidFunction(double min, double mid1, double mid2, double max);

  /**
   * Returns a function object that represents a Gaussian membership function.
   * @param mean The mean of the Gaussian. The membership function will be 1 for values equal to mean.
   * @param sigma The standard deviation of the Gaussian.
   */
  static std::function<double(double)> gaussianFunction(double mean, double sigma);

  /**
   * Returns a function object that represents a sigmoid membership function.
   * @param center The center of the sigmoid. The membership function will be 0.5 for values equal to center.
   * @param slope The slope of the sigmoid.
   */
  static std::function<double(double)> sigmoidFunction(double center, double slope);

  /**
   * Returns a function object that represents a finite sigmoid membership function.
   * @param lower The lower border of the sigmoid. The membership function will be 0 for values smaller than lower.
   * @param center The center of the sigmoid. The membership function will be 0.5 for values equal to center.
   * @param upper The upper border of the sigmoid. The membership function will be 1 for values greater than upper.
   */
  static std::function<double(double)> sigmoidFiniteFunction(double lower, double center, double upper);

  /**
   * Throws an exception with a message that the number of parameters is invalid.
   * @param linguisticTerm The linguistic term for which the FuzzySet should be created.
   * @param functionName The name of the function that should be created.
   * @param expected The expected number of parameters.
   * @param actual The actual number of parameters.
   */
  static void throwInvalidNumberOfArguments(const std::string &linguisticTerm, const std::string &functionName,
                                            size_t expected, size_t actual);
};

}  // namespace autopas::fuzzy_logic