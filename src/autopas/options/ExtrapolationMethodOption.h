/**
 * @file ExtrapolationMethodOption.h
 * @author Julian Pelloth
 * @date 06.05.2020
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the choices of possible extrapolation methods for predictiveTuning.
 */
class ExtrapolationMethodOption : public Option<ExtrapolationMethodOption> {
 public:
  /**
   * Possible choices for the auto tuner.
   */
  enum Value {
    /**
     * Places a line through the last two data points.
     */
    linePrediction,
    /**
     * Places a line through the last _evidenceFirstPrediction points.
     */
    linearRegression,
    /**
     * Places a polynomial function through a certain number of data points using Lagrange interpolation.
     */
    lagrange,
    /**
     * Places a polynomial function through a certain number of data points using Newtons method.
     */
    newton,
  };

  /**
   * Constructor.
   */
  ExtrapolationMethodOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr ExtrapolationMethodOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of ExtrapolationMethods.
   * @return map option -> string representation
   */
  static std::map<ExtrapolationMethodOption, std::string> getOptionNames() {
    return {
        {ExtrapolationMethodOption::linePrediction, "line-prediction"},
        {ExtrapolationMethodOption::linearRegression, "linear-regression"},
        {ExtrapolationMethodOption::lagrange, "lagrange"},
        {ExtrapolationMethodOption::newton, "newton"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
