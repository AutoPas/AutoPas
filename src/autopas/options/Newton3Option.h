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
   * Possible choices for the Newton3 optimization. Values consistent with bool.
   */
  enum Value {
    /**
     * Calculate F_{ij} and F_{ji} individually.
     */
    disabled,
    /**
     * Use F_{ij} = -F_{ji}.
     */
    enabled
  };

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
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<Newton3Option> getDiscouragedOptions() { return {}; }

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
