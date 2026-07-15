/**
 * @file SortingThresholdInfo2B.h
 * @date 14 Jul 2026
 */

#pragma once

#include <cstddef>

#include "autopas/options/CellLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SortingThresholdInfoInterface.h"

namespace autopas {

/**
 * Per-Newton3-state, per-CellLayoutOption pair-sorting thresholds for a 2-body CellFunctor.
 *
 * Replaces the previous opaque std::array<std::array<size_t, 3>, 2> grid with named members, accessed through
 * getThreshold()/setThreshold().
 */
class SortingThresholdInfo2B : public SortingThresholdInfoInterface {
 public:
  /**
   * Constructor. All thresholds are initialized to their in-class default (8).
   */
  SortingThresholdInfo2B() = default;

  /**
   * Constructor setting every Newton3-state / CellLayoutOption combination to the same value.
   * @param uniformThreshold Threshold value applied to all six combinations.
   */
  explicit SortingThresholdInfo2B(size_t uniformThreshold) {
    for (auto &row : _thresholds) {
      std::ranges::fill(row, uniformThreshold);
    }
  }

  /**
   * Set the threshold for a given Newton3 state and CellLayoutOption.
   * @param newton3
   * @param layout
   * @param value
   */
  void setThreshold(Newton3Option newton3, CellLayoutOption layout, size_t value) {
    _thresholds[static_cast<std::size_t>(newton3)][static_cast<std::size_t>(layout)] = value;
  }

  /**
   * Get the threshold for a given Newton3 state and CellLayoutOption.
   * @param newton3
   * @param layout
   * @return The stored threshold.
   */
  [[nodiscard]] size_t getThreshold(Newton3Option newton3, CellLayoutOption layout) const {
    return _thresholds[static_cast<std::size_t>(newton3)][static_cast<std::size_t>(layout)];
  }

  /**
   * Set the threshold for a given Newton3 state and CellLayoutOption.
   * @param newton3
   * @param layout
   * @param value
   */
  void setThreshold(bool newton3, std::array<double, 3> sortingDirection, size_t value) {
    setThreshold(newton3 ? Newton3Option::enabled : Newton3Option::disabled,
                 CellLayoutOption::fromSortingDirection(sortingDirection), value);
  }
  /**
   * Get the threshold for a given Newton3 state and CellLayoutOption.
   * @param newton3
   * @param layout
   * @return The stored threshold.
   */
  [[nodiscard]] size_t getThreshold(bool newton3, std::array<double, 3> sortingDirection) const {
    return getThreshold(newton3 ? Newton3Option::enabled : Newton3Option::disabled,
                        CellLayoutOption::fromSortingDirection(sortingDirection));
  }

 private:
  std::array<std::array<size_t, 3>, 2> _thresholds{};
};

}  // namespace autopas
