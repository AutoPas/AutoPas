/**
 * @file FuzzySetFactory.h
 * @author Manuel Lerchner
 * @date 19.04.24
 */

#pragma once

#include <functional>
#include <set>

#include "FuzzySet.h"

namespace autopas::FuzzyLogic {

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
   * Constructs a FuzzySet with the given linguistic term, based on the given activation function.
   * @param linguisticTerm The linguistic term of this FuzzySet.
   * @param functionName The name of the function to create.
   * @param params The parameters of the function.
   * @return A shared pointer to the created FuzzySet.
   */
  static std::shared_ptr<FuzzySet> makeFuzzySet(const std::string &linguisticTerm, const std::string &functionName,
                                                const std::vector<double> &params);

 private:
  /**
   * Returns a function object that represents a triangular membership function.
   * @param min The left border of the triangle. The membership function will be 0 for values smaller than min.
   * @param peak The peak of the triangle. The membership function will be 1 for values equal to peak.
   * @param max The right border of the triangle. The membership function will be 0 for values greater than max.
   */
  static std::function<double(double)> triangleFunction(double min, double peak, double max);

  /**
   * Returns a function object that represents a trapezoidal membership function.
   * @param min The left border of the trapezoid. The membership function will be 0 for values smaller than min.
   * @param leftPeak The left peak of the trapezoid. The membership function will be 1 for values in between leftPeak
   * and rightPeak.
   * @param rightPeak The right peak of the trapezoid. The membership function will be 1 for values in between leftPeak
   * and rightPeak.
   * @param max The right border of the trapezoid. The membership function will be 0 for values greater than max.
   */
  static std::function<double(double)> trapezoidFunction(double min, double leftPeak, double rightPeak, double max);

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
   *
   * The sigmoid function is defined as 1 / (1 + exp(-slope * (x - center))).
   *
   * The slope determines how steep the sigmoid is around the center. The derivative of the sigmoid at x=center is
   * given by slope / 4.
   */
  static std::function<double(double)> sigmoidFunction(double center, double slope);

  /**
   * Returns a function object that represents a finite sigmoid membership function.
   * @param lower The border of the sigmoid beyond which the membership function is 0. Does not necessarily mean the
   * border with the lowest x value.
   * @param center The center of the sigmoid. The membership function will be 0.5 for values equal to center.
   * @param upper The border of the sigmoid beyond which the membership function is 1. Does not necessarily mean the
   * border with the highest x value.
   *
   * Note: this is not an actual sigmoid function, but has an s-like shape. It is built using two parabolas and forms a
   * kind of cubic-shaped function.
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

/**
 * Class representing the choices for the membership functions.
 */
class MembershipFunctionOption : public Option<MembershipFunctionOption> {
 public:
  /**
   * Enum for the different defuzzification methods.
   */
  enum Value {
    /**
     * Triangle function
     */
    Triangle,
    /**
     * Trapezoid function
     */
    Trapezoid,
    /**
     * Gaussian function
     */
    Gaussian,
    /**
     * Sigmoid function
     */
    Sigmoid,
    /**
     * SigmoidFinite function
     */
    SigmoidFinite
  };

  /**
   * Constructor.
   */
  MembershipFunctionOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr MembershipFunctionOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of DefuzzificationMethod.
   * @return map option -> string representation
   */
  static std::map<MembershipFunctionOption, std::string> getOptionNames() {
    return {
        {MembershipFunctionOption::Triangle, "Triangle"},           {MembershipFunctionOption::Trapezoid, "Trapezoid"},
        {MembershipFunctionOption::Gaussian, "Gaussian"},           {MembershipFunctionOption::Sigmoid, "Sigmoid"},
        {MembershipFunctionOption::SigmoidFinite, "SigmoidFinite"},
    };
  };

 private:
  Value _value{Value(-1)};
};

}  // namespace autopas::FuzzyLogic