/**
 * @file SelectorStrategyOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <vector>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the choices for timing samples to be aggregated.
 */
class SelectorStrategyOption : public Option<SelectorStrategyOption> {
 public:
  /**
   * Possible choices for the employed selectors.
   */
  enum Value {
    /**
     * Fastest absolute value.
     */
    fastestAbs,
    /**
     * Fastest mean value.
     */
    fastestMean,
    /**
     * Fastest median value
     */
    fastestMedian
  };

  SelectorStrategyOption() = default;
  constexpr SelectorStrategyOption(Value option) : _value(option) {}
  constexpr operator Value() const { return _value; }
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   */
  static std::map<SelectorStrategyOption, std::string> getOptionNames() {
    return {
        {SelectorStrategyOption::fastestAbs, "Fastest-Absolute-Value"},
        {SelectorStrategyOption::fastestMean, "Fastest-Mean-Value"},
        {SelectorStrategyOption::fastestMedian, "Fastest-Median-Value"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas