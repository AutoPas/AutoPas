/**
 * @file TuningStrategyOption.h
 * @author F. Gratl
 * @date 6/3/19
 */

#pragma once

#include <set>

namespace autopas {

/**
 * Possible choices for the auto tuner.
 */
enum TuningStrategyOption {
  /**
   * Tests all allowed configurations and select the best.
   */
  exhaustiveSearch = 0,
};

/**
 * Provides a way to iterate over the possible choices of TuningStrategy.
 */
static const std::set<TuningStrategyOption> allTuningStrategyOptions = {
    TuningStrategyOption::exhaustiveSearch,
};
}  // namespace autopas