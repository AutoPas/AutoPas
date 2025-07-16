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
 * Class representing the different types of triggers that can be used to dynamically trigger tuning phases.
 */
class TuningTriggerOption : public Option<TuningTriggerOption> {
 public:
  /**
   * Possible choices for the auto tuner.
   */
  enum Value {
    /**
     * Default static trigger. Triggers tuning phase on multiples of tuningInterval.
     */
    staticSimple,
    /**
     * Compare runtime of current iteration to runtime of last iteration, retune if current iteration runtime is greater
     * or equal than triggerFactor times last iteration runtime.
     */
    timeBasedSimple,

    /**
     *  Similar to timeBasedSimple, but compares to the moving average of last n iterations.
     **/
    timeBasedAverage,
    /**
     *  Splits up last n samples into two intervals A, B. Triggers if avg(B) is greater or equal than triggerFactor times avg(A).
     **/
    timeBasedSplit,
    /**
     *  Performs a linear regression on the last n iteration runtimes. Triggers if the slope is bigger than triggerFactor.
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
