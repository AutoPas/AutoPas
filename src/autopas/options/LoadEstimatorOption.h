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
 * Class representing the traversal choices.
 */
class LoadEstimatorOption : public Option<LoadEstimatorOption> {
 public:
  /**
   * Possible choices for the cell pair traversal.
   */
  enum Value {
    none = 0,
    squaredParticlesPerCell = 1,
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
   * Provides a way to iterate over the possible choices of LoadEstimatorOption.
   * @return map option -> string representation
   */
  static std::map<LoadEstimatorOption, std::string> getOptionNames() {
    return {
        {LoadEstimatorOption::none, "none"},
        {LoadEstimatorOption::squaredParticlesPerCell, "squared-particles-per-cell"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
