/**
 * @file SelectorStrategie.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <vector>

namespace autopas {

/**
 * Possible choices for the employed selectors.
 */
enum SelectorStrategy {
  /**
   * Fastest absolute value.
   */
  fastestAbs,
  /**
   * Fastest mean value.
   */
  fastestMean,
  /**
   * Fastest median value
   */
  fastestMedian
};

/**
 * Provides a way to iterate over the possible choices of selector strategies.
 */
static const std::set<SelectorStrategy> allSelectorStrategies = {
    SelectorStrategy::fastestAbs,
    SelectorStrategy::fastestMean,
    SelectorStrategy::fastestMedian,
};
}  // namespace autopas