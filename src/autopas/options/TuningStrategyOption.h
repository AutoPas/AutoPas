/**
 * @file TuningStrategyOption.h
 * @author F. Gratl
 * @date 6/3/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the choices of possible tuning strategies for the auto-tuner.
 */
class TuningStrategyOption : public Option<TuningStrategyOption> {
 public:
  /**
   * Possible choices for the auto tuner.
   */
  enum Value {
    /**
     *  Random test configurations and select the best.
     **/
    randomSearch,
    /**
     * Tests all allowed configurations and select the best.
     */
    fullSearch,
    /**
     * Predict the configuration which will yield the most
     * information if tested next.
     */
    bayesianSearch,
    /**
     * ActiveHarmony client / server system
     */
    activeHarmony,

    /**
     * Predicts all configurations, tests those which are in the optimum range and selects the best.
     */
    predictiveTuning,
  };

  /**
   * Constructor.
   */
  TuningStrategyOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr TuningStrategyOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of TuningStrategy.
   * @return map option -> string representation
   */
  static std::map<TuningStrategyOption, std::string> getOptionNames() {
    return {
        {TuningStrategyOption::bayesianSearch, "bayesian-Search"},
        {TuningStrategyOption::fullSearch, "full-Search"},
        {TuningStrategyOption::randomSearch, "random-Search"},
        {TuningStrategyOption::activeHarmony, "active-harmony"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
