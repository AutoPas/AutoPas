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
   * Threshold for config {n3: off, cellDirection: Face, bidirectional: true}
   */
  size_t noN3FaceThresholdBidirectional;
  /**
   * Threshold for config {n3: off, cellDirection: Face, bidirectional: false}
   */
  size_t noN3FaceThresholdUnidirectional;
  /**
   * Threshold for config {n3: off, cellDirection: Edge, bidirectional: true}
   */
  size_t noN3EdgeThresholdBidirectional;
  /**
   * Threshold for config {n3: off, cellDirection: Edge, bidirectional: false}
   */
  size_t noN3EdgeThresholdUnidirectional;
  /**
   * Threshold for config {n3: off, cellDirection: Corner, bidirectional: true}
   */
  size_t noN3CornerThresholdBidirectional;
  /**
   * Threshold for config {n3: off, cellDirection: Corner, bidirectional: false}
   */
  size_t noN3CornerThresholdUnidirectional;
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
   * Constructor setting every Newton3-state / SortingDirectionOption / bidirectional combination to the same
   * value.
   * @param uniformThreshold Threshold value applied to all nine combinations.
   */
  explicit SortingThresholdInfo2B(size_t uniformThreshold)
      : noN3FaceThresholdBidirectional(uniformThreshold),
        noN3FaceThresholdUnidirectional(uniformThreshold),
        noN3EdgeThresholdBidirectional(uniformThreshold),
        noN3EdgeThresholdUnidirectional(uniformThreshold),
        noN3CornerThresholdBidirectional(uniformThreshold),
        noN3CornerThresholdUnidirectional(uniformThreshold),
        n3FaceThreshold(uniformThreshold),
        n3EdgeThreshold(uniformThreshold),
        n3CornerThreshold(uniformThreshold) {}

  /**
   * Constructor to set each per Newton3/CellDirection/bidirectional combination value explicitly. Newton3-on
   * values are not split by bidirectional, since bidirectional only affects the measured timing when Newton3
   * is off.
   * @param noN3FaceThresholdBidirectional
   * @param noN3FaceThresholdUnidirectional
   * @param noN3EdgeThresholdBidirectional
   * @param noN3EdgeThresholdUnidirectional
   * @param noN3CornerThresholdBidirectional
   * @param noN3CornerThresholdUnidirectional
   * @param n3FaceThreshold
   * @param n3EdgeThreshold
   * @param n3CornerThreshold
   */
  SortingThresholdInfo2B(size_t noN3FaceThresholdBidirectional, size_t noN3FaceThresholdUnidirectional,
                         size_t noN3EdgeThresholdBidirectional, size_t noN3EdgeThresholdUnidirectional,
                         size_t noN3CornerThresholdBidirectional, size_t noN3CornerThresholdUnidirectional,
                         size_t n3FaceThreshold, size_t n3EdgeThreshold, size_t n3CornerThreshold)
      : noN3FaceThresholdBidirectional(noN3FaceThresholdBidirectional),
        noN3FaceThresholdUnidirectional(noN3FaceThresholdUnidirectional),
        noN3EdgeThresholdBidirectional(noN3EdgeThresholdBidirectional),
        noN3EdgeThresholdUnidirectional(noN3EdgeThresholdUnidirectional),
        noN3CornerThresholdBidirectional(noN3CornerThresholdBidirectional),
        noN3CornerThresholdUnidirectional(noN3CornerThresholdUnidirectional),
        n3FaceThreshold(n3FaceThreshold),
        n3EdgeThreshold(n3EdgeThreshold),
        n3CornerThreshold(n3CornerThreshold) {}

  /**
   * Getter to automatically get the correct threshold based on current newton3/CellDirection/bidirectional
   * configuration. bidirectional is ignored when newton3 is enabled, since it does not affect the timing of
   * the newton3-on path (see SortingThresholdBenchmark).
   * @param newton3 whether the newton3 optimization is enabled or not.
   * @param sortingDirection the cell axis used for sorting. CellDirection is deduced from this.
   * @param bidirectional whether the CellFunctor processing this cell pair needs to compute both interaction
   * directions itself. Only relevant when newton3 is disabled.
   * @return correct threshold value for the given configuration.
   */
  size_t getThresholdByConfig(bool newton3, std::array<double, 3> sortingDirection, bool bidirectional) const {
    const SortingDirectionOption cellDirection = SortingDirectionOption::fromRawArray(sortingDirection);
    if (newton3) {
      if (cellDirection == SortingDirectionOption::corner) {
        return n3CornerThreshold;
      }
      if (cellDirection == SortingDirectionOption::edge) {
        return n3EdgeThreshold;
      }
      return n3FaceThreshold;
    } else if (bidirectional) {
      if (cellDirection == SortingDirectionOption::corner) {
        return noN3CornerThresholdBidirectional;
      }
      if (cellDirection == SortingDirectionOption::edge) {
        return noN3EdgeThresholdBidirectional;
      }
      return noN3FaceThresholdBidirectional;
    } else {
      if (cellDirection == SortingDirectionOption::corner) {
        return noN3CornerThresholdUnidirectional;
      }
      if (cellDirection == SortingDirectionOption::edge) {
        return noN3EdgeThresholdUnidirectional;
      }
      return noN3FaceThresholdUnidirectional;
    }
  }
  /**
   * Setter to set struct values by configuration, represented by options. bidirectional is ignored when n3 is
   * Newton3Option::enabled, since it does not affect the timing of the newton3-on path (see
   * SortingThresholdBenchmark).
   * @param n3 Whether the Newton3 optimization is used.
   * @param cellDirection How the cells are placed in relation to each other.
   * @param bidirectional Whether the CellFunctor processing this cell pair needs to compute both interaction
   * directions itself. Only relevant when n3 is Newton3Option::disabled.
   * @param value The value to set.
   */
  void setThresholdByOption(Newton3Option n3, SortingDirectionOption cellDirection, bool bidirectional, size_t value) {
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
    } else if (bidirectional) {
      switch (cellDirection) {
        case SortingDirectionOption::corner:
          noN3CornerThresholdBidirectional = value;
          break;
        case SortingDirectionOption::edge:
          noN3EdgeThresholdBidirectional = value;
          break;
        case SortingDirectionOption::face:
          noN3FaceThresholdBidirectional = value;
          break;
      }
    } else {
      switch (cellDirection) {
        case SortingDirectionOption::corner:
          noN3CornerThresholdUnidirectional = value;
          break;
        case SortingDirectionOption::edge:
          noN3EdgeThresholdUnidirectional = value;
          break;
        case SortingDirectionOption::face:
          noN3FaceThresholdUnidirectional = value;
          break;
      }
    }
  }
};

}  // namespace autopas
