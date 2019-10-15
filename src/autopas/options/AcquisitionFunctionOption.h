/**
 * @file AcquisitionFunctionOption.h
 * @author Jan Nguyen
 * @date 26.06.19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the acquisition function choices for the Bayesian search.
 */
class AcquisitionFunctionOption : public Option<AcquisitionFunctionOption> {
 public:
  /**
   * Different acquisition functions
   */
  enum Value {
    /**
     * Upper confidence bound
     */
    ucb,
    /**
     * Lower confidence bound
     */
    lcb,
    /**
     * mean
     */
    mean
  };

  AcquisitionFunctionOption() = default;
  constexpr AcquisitionFunctionOption(Value option) : _value(option) {}
  constexpr operator Value() const { return _value; }
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   */
  static std::map<AcquisitionFunctionOption, std::string> getOptionNames() {
    return {
        {AcquisitionFunctionOption::ucb, "upper-confidence-bound"},
        {AcquisitionFunctionOption::lcb, "lower-confidence-bound"},
        {AcquisitionFunctionOption::mean, "mean"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
