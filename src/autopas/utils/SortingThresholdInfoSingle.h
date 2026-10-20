/**
 * @file SortingThresholdInfoSingle.h
 * @date 14.07.2026
 * @author hmeyran
 */

#pragma once

#include <cstddef>

#include "autopas/utils/SortingThresholdInfoInterface.h"

namespace autopas {

/**
 * Single threshold value struct.
 */
struct SortingThresholdInfoSingle : SortingThresholdInfoInterface {
  /**
   * Constructor for single threshold value.
   * @param threshold
   */
  SortingThresholdInfoSingle(size_t threshold) : threshold(threshold) {}
  /**
   * Actual threshold data.
   */
  size_t threshold;
};

}  // namespace autopas
