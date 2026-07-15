/**
 * @file SortingThresholdInfoSingle.h
 * @date 14 Jul 2026
 */

#pragma once

#include <cstddef>

#include "autopas/utils/SortingThresholdInfoInterface.h"

namespace autopas {

/**
 * Per-Newton3-state, per-CellLayoutOption pair-sorting thresholds for a 2-body CellFunctor.
 *
 * Replaces the previous opaque std::array<std::array<size_t, 3>, 2> grid with named members, accessed through
 * getThreshold()/setThreshold().
 */
class SortingThresholdInfoSingle : public SortingThresholdInfoInterface {
 public:
  /**
   * Constructor. All thresholds are initialized to their in-class default (8).
   */
  SortingThresholdInfoSingle() = default;

  /**
   * Constructor setting every Newton3-state / CellLayoutOption combination to the same value.
   * @param threshold Threshold value applied to all six combinations.
   */
  explicit SortingThresholdInfoSingle(size_t threshold) : _threshold(threshold) {}

  /**
   * Set the threshold for a given Newton3 state and CellLayoutOption.
   * @param value
   */
  void setThreshold(size_t value) { _threshold = value; }
  /**
   * Get the threshold for a given Newton3 state and CellLayoutOption.

   * @return The stored threshold.
   */
  [[nodiscard]] size_t getThreshold() const { return _threshold; }

 private:
  size_t _threshold;
};

}  // namespace autopas
