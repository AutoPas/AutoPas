/**
 * @file DistanceClassOption.h
 * @author D. Martin
 * @date 05/06/2025
 */

#pragma once

#include <set>

#include "Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the choices for the Newton 3 optimization.
 */
class DistanceClassOption : public Option<DistanceClassOption> {
 public:
  /**
   * Possible choices for the Newton3 optimization. Values consistent with bool.
   */
  enum Value {
    /**
     * Use IBI force in outer distance class.
     */
    ibi,
    /**
     * Use full particle force in outer distance class.
     */
    fp,
    /**
     * Use coarse grain particle force in outer distance class.
     */
    cgmol,
    /**
     * Disable distance class
     */
    disabled,
  };

  /**
   * Constructor.
   */
  DistanceClassOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr DistanceClassOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<DistanceClassOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   * @return map option -> string representation
   */
  static std::map<DistanceClassOption, std::string> getOptionNames() {
    return {
        {DistanceClassOption::ibi, "ibi"},
        {DistanceClassOption::fp, "fp"},
        {DistanceClassOption::cgmol, "cgmol"},
        {DistanceClassOption::disabled, "disabled"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
