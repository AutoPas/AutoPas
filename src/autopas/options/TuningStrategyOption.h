/**
 * @file TuningStrategyOption.h
 * @author F. Gratl
 * @date 03.06.2019
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
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
     * Technically this is no tuning strategy anymore, but the option still exists for compatibility reasons.
     */
    fullSearch,
    /**
     * Predict the configuration which will yield the most
     * information if tested next.
     */
    bayesianSearch,
    /**
     * Predict the configuration which will yield the most
     * information if tested next. Uses a Gaussian Process
     * per allowed discrete tuple.
     */
    bayesianClusterSearch,
    /**
     * ActiveHarmony client / server system
     */
    activeHarmony,
    /**
     * Predicts performance of all configurations based on previous tuning phases, tests those which are in the optimum
     * range, and selects the best.
     */
    predictiveTuning,
    /**
     * Applies predefined rules to dynamically exclude configurations from tuning that are expected to perform worse
     * than others in the next tuning phase.
     */
    ruleBasedTuning,
    /**
     * Dynamic blacklist that throws out configurations that perform poorly.
     */
    slowConfigFilter,
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
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<TuningStrategyOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of TuningStrategy.
   * @return map option -> string representation
   */
  static std::map<TuningStrategyOption, std::string> getOptionNames() {
    return {
        {TuningStrategyOption::bayesianSearch, "bayesian-Search"},
        {TuningStrategyOption::bayesianClusterSearch, "bayesian-cluster-Search"},
        {TuningStrategyOption::fullSearch, "full-Search"},
        {TuningStrategyOption::randomSearch, "random-Search"},
        {TuningStrategyOption::activeHarmony, "active-harmony"},
        {TuningStrategyOption::predictiveTuning, "predictive-tuning"},
        {TuningStrategyOption::ruleBasedTuning, "rule-based-tuning"},
        {TuningStrategyOption::ruleBasedTuning, "slow-config-filter"},
    };
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
