/**
 * @file SelectorStrategies.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

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

}  // namespace autopas