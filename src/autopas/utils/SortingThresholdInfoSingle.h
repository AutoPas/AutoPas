/**
 * @file SortingThresholdInfoSingle.h
 * @date 14 Jul 2026
 * @author hmeyran
 */

#pragma once

#include <cstddef>

#include "autopas/utils/SortingThresholdInfoInterface.h"

namespace autopas {

/**
 * Single threshold value struct.
 * Stand in for the old single value way of handling the sorting Threshold. Used as a fallback option and the default
 * value.
 */
struct SortingThresholdInfoSingle : SortingThresholdInfoInterface {
  /**
   * Constructor for single threshold value.
   * Intentionally not explicit to make it easier to use this like you would the threshold before.
   * @param threshold
   */
  SortingThresholdInfoSingle(size_t threshold) : threshold(threshold) {}
  /**
   * Actual threshold data.
   */
  size_t threshold;
};

}  // namespace autopas
