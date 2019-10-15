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
     * Tests all allowed configurations and select the best.
     */
    fullSearch = 0,
    /**
     * Predict the configuration which will yield the most
     * information if tested next.
     */
    bayesianSearch = 1
  };

  TuningStrategyOption() = default;
  constexpr TuningStrategyOption(Value option) : _value(option) {}
  constexpr operator Value() const { return _value; }
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   */
  static std::map<TuningStrategyOption, std::string> getOptionNames() {
    return {
        {TuningStrategyOption::fullSearch, "full-Search"},
        {TuningStrategyOption::bayesianSearch, "bayesian-Search"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
