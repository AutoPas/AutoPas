/**
 * @file OptimumSelector.h
 * @author F. Gratl
 * @date 3/15/19
 */

#pragma once

#include <algorithm>
#include <numeric>
#include <vector>
#include "autopas/options/SelectorStrategie.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Collection of functions for selecting the optimum value out of a vector of values according to a given strategy
 */
namespace OptimumSelector {

/**
 * Minimal value.
 * @param values
 * @return Smallest value of the vector.
 */
inline unsigned long minValue(std::vector<unsigned long> values) {
  return *std::min_element(values.begin(), values.end());
}

/**
 * Arithmetic mean.
 * @param values
 * @return Arithmetic mean of the vector.
 */
inline unsigned long meanValue(std::vector<unsigned long> values) {
  return std::accumulate(values.begin(), values.end(), 0l) / values.size();
}

/**
 * Median value.
 * @param values
 * @return Middle ((size-1) /2) of the sorted vector.
 */
inline unsigned long medianValue(std::vector<unsigned long> values) {
  if (values.empty()) return 0;

  std::sort(values.begin(), values.end());
  return values[(values.size() - 1) / 2];
}

/**
 * Optimal value according to passed strategy.
 * @param values
 * @param strategy For possible selector strategy choices see AutoPas::SelectorStrategy.
 * @return value
 */
inline unsigned long optimumValue(const std::vector<unsigned long>& values, SelectorStrategy strategy) {
  switch (strategy) {
    case SelectorStrategy::fastestAbs: {
      return minValue(values);
    }
    case SelectorStrategy::fastestMean: {
      return meanValue(values);
    }
    case SelectorStrategy::fastestMedian: {
      return medianValue(values);
    }
      // no default to get appropriate compiler warnings
  }

  autopas::utils::ExceptionHandler::exception("OptionSelector: Unknown selector strategy {}!", strategy);

  return 0;
}

}  // namespace OptimumSelector
}  // namespace autopas
