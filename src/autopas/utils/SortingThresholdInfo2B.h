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
struct SortingThresholdInfo2B : SortingThresholdInfoInterface {
  size_t noN3FaceThreshold;
  size_t noN3EdgeThreshold;
  size_t noN3CornerThreshold;

  size_t n3FaceThreshold;
  size_t n3EdgeThreshold;
  size_t n3CornerThreshold;

  /**
   * Constructor setting every Newton3-state / CellLayoutOption combination to the same value.
   * @param uniformThreshold Threshold value applied to all six combinations.
   */
  explicit SortingThresholdInfo2B(size_t uniformThreshold)
      : n3FaceThreshold(uniformThreshold),
        n3EdgeThreshold(uniformThreshold),
        n3CornerThreshold(uniformThreshold),
        noN3FaceThreshold(uniformThreshold),
        noN3EdgeThreshold(uniformThreshold),
        noN3CornerThreshold(uniformThreshold) {}

  SortingThresholdInfo2B(size_t noN3FaceThreshold, size_t noN3EdgeThreshold, size_t noN3CornerThreshold,
                         size_t n3FaceThreshold, size_t n3EdgeThreshold, size_t n3CornerThreshold)
      : noN3FaceThreshold(noN3FaceThreshold),
        noN3EdgeThreshold(noN3EdgeThreshold),
        noN3CornerThreshold(noN3CornerThreshold),
        n3FaceThreshold(n3FaceThreshold),
        n3EdgeThreshold(n3EdgeThreshold),
        n3CornerThreshold(n3CornerThreshold) {}

  size_t getThresholdByConfig(bool newton3, std::array<double, 3> sortingDirection) const {
    size_t zeroes = std::ranges::count(sortingDirection, 0);
    if (newton3) {
      if (zeroes == 0) {
        return n3CornerThreshold;
      }
      if (zeroes == 1) {
        return n3EdgeThreshold;
      }
      return n3FaceThreshold;
    } else {
      if (zeroes == 0) {
        return noN3CornerThreshold;
      }
      if (zeroes == 1) {
        return noN3EdgeThreshold;
      }
      return noN3FaceThreshold;
    }
  }
};

}  // namespace autopas
