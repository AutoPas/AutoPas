/**
 * @file Newton3Option.h
 * @author F. Gratl
 * @date 01/02/2019
 */

#pragma once

#include <set>

#include "Option.h"

namespace autopas {

/**
 * Class representing the choices for the Newton 3 optimization.
 */
class Newton3Option : public Option<Newton3Option> {
 public:
  /**
   * Possible choices for the particle data layout. Values consistent with bool.
   */
  enum Value { disabled = 0, enabled = 1 };

  Newton3Option() = default;
  constexpr Newton3Option(Value option) : _value(option) {}
  constexpr operator Value() const { return _value; }
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
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
}  // namespace autopas
