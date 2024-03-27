/**
 * @file TuningMetricOption.h
 * @author fischerv
 * @date 09 May 2020
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the load estimator choices.
 */
class TuningMetricOption : public Option<TuningMetricOption> {
 public:
  /**
   * Possible choices for the load estimation algorithm.
   */
  enum Value {
    /**
     * Optimize for shortest simulation time
     */
    time,
    /**
     * Optimize for least energy usage
     */
    energy,
  };

  /**
   * Constructor.
   */
  TuningMetricOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr TuningMetricOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<TuningMetricOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of LoadEstimatorOption.
   * @return map option -> string representation
   */
  static std::map<TuningMetricOption, std::string> getOptionNames() {
    return {
        {TuningMetricOption::time, "time"},
        {TuningMetricOption::energy, "energy"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
