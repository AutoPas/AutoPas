/**
 * @file SortingThresholdInfo2B.h
 * @date 14.07.2026
 * @author hmeyran
 */

#pragma once

#include <cstddef>

#include "autopas/options/CellLayoutOption.h"
#include "autopas/utils/SortingThresholdInfoInterface.h"

namespace autopas {

/**
 * Per-Newton3-state, per-CellLayoutOption pair-sorting thresholds for a 2-body CellFunctor.
 * Provides publicly accessible named data as well as a getter to easily get the correct Value from a configuration.
 */
struct SortingThresholdInfo2B : SortingThresholdInfoInterface {
  /**
   * Threshold for config {n3: off, cellLayout: Face}
   */
  size_t noN3FaceThreshold;
  /**
   * Threshold for config {n3: off, cellLayout: Edge}
   */
  size_t noN3EdgeThreshold;
  /**
   * Threshold for config {n3: off, cellLayout: Corner}
   */
  size_t noN3CornerThreshold;
  /**
   * Threshold for config {n3: on, cellLayout: Face}
   */
  size_t n3FaceThreshold;
  /**
   * Threshold for config {n3: on, cellLayout: Edge}
   */
  size_t n3EdgeThreshold;
  /**
   * Threshold for config {n3: on, cellLayout: Corner}
   */
  size_t n3CornerThreshold;

  /**
   * Constructor setting every Newton3-state / CellLayoutOption combination to the same value.
   * @param uniformThreshold Threshold value applied to all six combinations.
   */
  explicit SortingThresholdInfo2B(size_t uniformThreshold)
      : noN3FaceThreshold(uniformThreshold),
        noN3EdgeThreshold(uniformThreshold),
        noN3CornerThreshold(uniformThreshold),
        n3FaceThreshold(uniformThreshold),
        n3EdgeThreshold(uniformThreshold),
        n3CornerThreshold(uniformThreshold) {}
  /**
   * Constructor to set each per Newton3/CellLayout combination value explicitly.
   * @param noN3FaceThreshold
   * @param noN3EdgeThreshold
   * @param noN3CornerThreshold
   * @param n3FaceThreshold
   * @param n3EdgeThreshold
   * @param n3CornerThreshold
   */
  SortingThresholdInfo2B(size_t noN3FaceThreshold, size_t noN3EdgeThreshold, size_t noN3CornerThreshold,
                         size_t n3FaceThreshold, size_t n3EdgeThreshold, size_t n3CornerThreshold)
      : noN3FaceThreshold(noN3FaceThreshold),
        noN3EdgeThreshold(noN3EdgeThreshold),
        noN3CornerThreshold(noN3CornerThreshold),
        n3FaceThreshold(n3FaceThreshold),
        n3EdgeThreshold(n3EdgeThreshold),
        n3CornerThreshold(n3CornerThreshold) {}
  /**
   * Getter to automatically get the correct threshold based on current newton3/CellLayout configuration
   * @param newton3 whether the newton3 optimization is enabled or not.
   * @param sortingDirection the cell axis used for sorting. CellLayout is deduced from this.
   * @return correct threshold value for the given configuration.
   */
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
