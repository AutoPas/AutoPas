/**
 * @file CellLayoutOption.h
 * @date 14.07.2026
 * @author hmeyran
 */

#pragma once

#include <algorithm>
#include <array>
#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the classification of a normalized cell-pair sortingDirection by the number of its zero
 * components.
 *
 * // TODO: name collides conceptually with DataLayoutOption (AoS/SoA); rename if a clearer name is found
 * // (candidates considered: SortingDirectionOption, DirectionTypeOption).
 */
class CellLayoutOption : public Option<CellLayoutOption> {
 public:
  /**
   * Possible choices for the direction-type classification, indexed by the number of zero components in the
   * normalized sortingDirection.
   */
  enum Value : size_t {
    /**
     * All three components of sortingDirection are non-zero.
     */
    corner = 0,
    /**
     * Exactly one component of sortingDirection is zero.
     */
    edge = 1,
    /**
     * Exactly two components of sortingDirection are zero.
     */
    face = 2
  };

  /**
   * Constructor.
   */
  CellLayoutOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr CellLayoutOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<CellLayoutOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of CellLayoutOption.
   * @return map option -> string representation
   */
  static std::map<CellLayoutOption, std::string> getOptionNames() {
    return {
        {CellLayoutOption::corner, "Corner"},
        {CellLayoutOption::edge, "Edge"},
        {CellLayoutOption::face, "Face"},
    };
  };

  /**
   * Classifies a normalized sortingDirection by its number of zero components.
   * @param sortingDirection Normalized vector connecting the centers of two cells. Must not be {0., 0., 0.}.
   * @return The corresponding CellLayoutOption.
   */
  static CellLayoutOption fromSortingDirection(const std::array<double, 3> &sortingDirection) {
    return static_cast<Value>(std::ranges::count(sortingDirection, 0.0));
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
