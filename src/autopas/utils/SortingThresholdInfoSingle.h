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
 * Per-Newton3-state, per-CellLayoutOption pair-sorting thresholds for a 2-body CellFunctor.
 *
 * Replaces the previous opaque std::array<std::array<size_t, 3>, 2> grid with named members, accessed through
 * getThreshold()/setThreshold().
 */
struct SortingThresholdInfoSingle : SortingThresholdInfoInterface {
  SortingThresholdInfoSingle(size_t threshold) : threshold(threshold) {}
  size_t threshold;
};

}  // namespace autopas
