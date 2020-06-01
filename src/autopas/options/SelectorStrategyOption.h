/**
 * @file SelectorStrategyOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <vector>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
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

  /**
   * Constructor.
   */
  SelectorStrategyOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr SelectorStrategyOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<SelectorStrategyOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   * @return map option -> string representation
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
}  // namespace options
}  // namespace autopas