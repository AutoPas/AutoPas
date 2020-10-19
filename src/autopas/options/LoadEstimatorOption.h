/**
 * @file LoadEstimatorOption.h
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
class LoadEstimatorOption : public Option<LoadEstimatorOption> {
 public:
  /**
   * Possible choices for the load estimation algorithm.
   */
  enum Value {
    /**
     * No load estimator. If the Configuration supports load estimators, everything is assigned the same estimation.
     */
    none,
    /**
     * Number of particles per cell squared.
     */
    squaredParticlesPerCell,
    /**
     * Sum of neighbor list lengths.
     */
    neighborListLength,
  };

  /**
   * Constructor.
   */
  LoadEstimatorOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr LoadEstimatorOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<LoadEstimatorOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of LoadEstimatorOption.
   * @return map option -> string representation
   */
  static std::map<LoadEstimatorOption, std::string> getOptionNames() {
    return {
        {LoadEstimatorOption::none, "none"},
        {LoadEstimatorOption::squaredParticlesPerCell, "squared-particles-per-cell"},
        {LoadEstimatorOption::neighborListLength, "neighbor-list-length"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
