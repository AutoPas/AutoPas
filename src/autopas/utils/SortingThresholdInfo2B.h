/**
 * @file SortingThresholdInfo2B.h
 * @date 14.07.2026
 * @author hmeyran
 */

#pragma once

#include <cstddef>

#include "autopas/options/Newton3Option.h"
#include "autopas/options/SortingDirectionOption.h"
#include "autopas/utils/SortingThresholdInfoInterface.h"

namespace autopas {

/**
 * Per-Newton3-state, per-SortingDirectionOption pair-sorting thresholds for a 2-body CellFunctor.
 * Provides publicly accessible named data as well as a getter to easily get the correct Value from a configuration.
 */
struct SortingThresholdInfo2B : SortingThresholdInfoInterface {
  /**
   * Threshold for config {n3: off, cellDirection: Face}
   */
  size_t noN3FaceThreshold;
  /**
   * Threshold for config {n3: off, cellDirection: Edge}
   */
  size_t noN3EdgeThreshold;
  /**
   * Threshold for config {n3: off, cellDirection: Corner}
   */
  size_t noN3CornerThreshold;
  /**
   * Threshold for config {n3: on, cellDirection: Face}
   */
  size_t n3FaceThreshold;
  /**
   * Threshold for config {n3: on, cellDirection: Edge}
   */
  size_t n3EdgeThreshold;
  /**
   * Threshold for config {n3: on, cellDirection: Corner}
   */
  size_t n3CornerThreshold;

  /**
   * Constructor setting every Newton3-state / SortingDirectionOption combination to the same value.
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
   * Constructor to set each per Newton3/CellDirection combination value explicitly.
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
   * Getter to automatically get the correct threshold based on current newton3/CellDirection configuration
   * @param newton3 whether the newton3 optimization is enabled or not.
   * @param sortingDirection the cell axis used for sorting. CellDirection is deduced from this.
   * @return correct threshold value for the given configuration.
   */
  size_t getThresholdByConfig(bool newton3, std::array<double, 3> sortingDirection) const {
    const SortingDirectionOption cellDirection = SortingDirectionOption::fromRawArray(sortingDirection);
    if (newton3) {
      if (cellDirection == SortingDirectionOption::corner) {
        return n3CornerThreshold;
      }
      if (cellDirection == SortingDirectionOption::edge) {
        return n3EdgeThreshold;
      }
      return n3FaceThreshold;
    } else {
      if (cellDirection == SortingDirectionOption::corner) {
        return noN3CornerThreshold;
      }
      if (cellDirection == SortingDirectionOption::edge) {
        return noN3EdgeThreshold;
      }
      return noN3FaceThreshold;
    }
  }
  /**
   * Setter to set struct values by configuration, represented by options.
   * @param n3 Whether the Newton3 optimization is used.
   * @param cellDirection How the cells are placed in relation to each other.
   * @param value The value to set.
   */
  void setThresholdByOption(Newton3Option n3, SortingDirectionOption cellDirection, size_t value) {
    if (n3 == Newton3Option::enabled) {
      switch (cellDirection) {
        case SortingDirectionOption::corner:
          n3CornerThreshold = value;
          break;
        case SortingDirectionOption::edge:
          n3EdgeThreshold = value;
          break;
        case SortingDirectionOption::face:
          n3FaceThreshold = value;
          break;
      }
    } else {
      switch (cellDirection) {
        case SortingDirectionOption::corner:
          noN3CornerThreshold = value;
          break;
        case SortingDirectionOption::edge:
          noN3EdgeThreshold = value;
          break;
        case SortingDirectionOption::face:
          noN3FaceThreshold = value;
          break;
      }
    }
  }
};

}  // namespace autopas
