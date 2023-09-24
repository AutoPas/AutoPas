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
class BoundaryTypeOption : public autopas::Option<BoundaryTypeOption> {
 public:
  /**
   * Possible choices of boundary condition
   */
  enum Value {
    /**
     * Periodic.
     */
    periodic,
    /**
     * Reflective.
     */
    reflective,
    /**
     * No Boundary.
     */
    none
  };
  BoundaryTypeOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr BoundaryTypeOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<BoundaryTypeOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of TuningStrategy.
   * @return map option -> string representation
   */
  static std::map<BoundaryTypeOption, std::string> getOptionNames() {
    return {{BoundaryTypeOption::periodic, "periodic"},
            {BoundaryTypeOption::reflective, "reflective"},
            {BoundaryTypeOption::none, "none"}};
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options