/**
 * @file DataLayoutOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the traversal choices.
 */
class DataLayoutOption : public Option<DataLayoutOption> {
 public:
  /**
   * Possible choices for the particle data layout.
   */
  enum Value {
    /**
     * Array of structure (AoS) : Particles are stored as Objects in a vector.
     */
    aos,
    /**
     * Structure of arrays (SoA) : One Vector for each particle property.
     */
    soa,
  };

  /**
   * Constructor.
   */
  DataLayoutOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr DataLayoutOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<DataLayoutOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   * @return map option -> string representation
   */
  static std::map<DataLayoutOption, std::string> getOptionNames() {
    return {
        {DataLayoutOption::aos, "AoS"},
        {DataLayoutOption::soa, "SoA"},
    };
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
