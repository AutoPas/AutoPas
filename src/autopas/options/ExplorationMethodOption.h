/**
 * @file ExplorationMethodOption.h
 * @author P. Metscher
 * @date 02.01.2026
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Option for the exploration method used in the Deep Reinforcement Learning tuning strategy.
 */
class ExplorationMethodOption : public Option<ExplorationMethodOption> {
 public:
  /**
   * Possible exploration methods for the Deep Reinforcement Learning tuning strategy.
   */
  enum Value {
    /**
     * The polynomial exploration method.
     *
     * Choose the next configuration based on the weighted age squared and the predicted evidence:
     *
     * @f[age^2 \cdot _phaseScale + predictedEvidence()@f]
     *
     * The predicted evidence is an interpolation of the last two tuning phases.
     */
    polynomial,
    /**
     * The random exploration method. Choose the next exploration samples randomly.
     */
    random,
    /**
     * The longest ago exploration method. Explore the configurations that have not been searched for the longest time.
     */
    longestAgo,
  };

  /**
   * Constructor.
   */
  ExplorationMethodOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr ExplorationMethodOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of ExplorationMethods.
   * @return map option -> string representation
   */
  static std::map<ExplorationMethodOption, std::string> getOptionNames() {
    return {
        {Value::polynomial, "polynomial"},
        {Value::random, "random"},
        {Value::longestAgo, "longest-ago"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas