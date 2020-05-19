/**
 * @file Newton3Option.h
 * @author F. Gratl
 * @date 01/02/2019
 */

#pragma once

#include <set>

#include "Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the choices for the Newton 3 optimization.
 */
class Newton3Option : public Option<Newton3Option> {
 public:
  /**
   * Possible choices for the particle data layout. Values consistent with bool.
   */
  enum Value { disabled = 0, enabled = 1 };

  /**
   * Constructor.
   */
  Newton3Option() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr Newton3Option(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   * @return map option -> string representation
   */
  static std::map<Newton3Option, std::string> getOptionNames() {
    return {
        {Newton3Option::disabled, "disabled"},
        {Newton3Option::enabled, "enabled"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
