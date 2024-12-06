/**
 * @file BoundaryTypeOption.h
 * @author S. J. Newcome
 * @date 24/01/2022
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace options {
/**
 * Class representing the options of boundary types
 */
class ZonalMethodOption : public autopas::Option<ZonalMethodOption> {
 public:
  /**
   * Possible choices of boundary condition
   */
  enum Value {
    fullshell,
    halfshell,
    midpoint
  };
  ZonalMethodOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr ZonalMethodOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<ZonalMethodOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of TuningStrategy.
   * @return map option -> string representation
   */
  static std::map<ZonalMethodOption, std::string> getOptionNames() {
    return {{ZonalMethodOption::fullshell, "fullshell"},
            {ZonalMethodOption::halfshell, "halfshell"},
            {ZonalMethodOption::midpoint, "midpoint"}};
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
