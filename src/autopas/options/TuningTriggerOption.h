/**
 * @file TuningTriggerOption.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the different types of triggers that can be used to dynamically initiate tuning phases.
 */
class TuningTriggerOption : public Option<TuningTriggerOption> {
 public:
  /**
   * Possible choices for the auto tuner.
   */
  enum Value {
    /**
     * Same behavior as with the default static tuning intervals. Triggers tuning phases on multiples of tuningInterval.
     */
    staticSimple,
    /**
     * Compares runtime of current iteration to runtime of last iteration, retunes if current iteration runtime is
     * greater than or equal to triggerFactor times last iteration runtime.
     */
    timeBasedSimple,

    /**
     * Similar to timeBasedSimple, but compares the current iteration runtime to the moving average of last n
     * iterations.
     **/
    timeBasedAverage,
    /**
     * Splits up last n runtime samples samples into two intervals A, B. Triggers if avg(B) is greater than or equal to
     * triggerFactor times avg(A).
     **/
    timeBasedSplit,
    /**
     * Performs a linear regression on the last n runtime samples. Triggers if the normalized slope estimator is
     * greater than or equal to triggerFactor.
     **/
    timeBasedRegression,
  };

  /**
   * Constructor.
   */
  TuningTriggerOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr TuningTriggerOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of TuningStrategy.
   * @return map option -> string representation
   */
  static std::map<TuningTriggerOption, std::string> getOptionNames() {
    return {
        {TuningTriggerOption::staticSimple, "StaticSimple"},
        {TuningTriggerOption::timeBasedSimple, "TimeBasedSimple"},
        {TuningTriggerOption::timeBasedAverage, "TimeBasedAverage"},
        {TuningTriggerOption::timeBasedSplit, "TimeBasedSplit"},
        {TuningTriggerOption::timeBasedRegression, "TimeBasedRegression"},
    };
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
