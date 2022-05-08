/**
 * @file OptimumSelector.h
 * @author F. Gratl
 * @date 3/15/19
 */

#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/utils/ExceptionHandler.h"

/**
 * Collection of functions for selecting the optimum value out of a vector of values according to a given strategy
 */
namespace autopas::OptimumSelector {

/**
 * Minimal value.
 * @param values
 * @return Smallest value of the vector.
 */
inline long minValue(const std::vector<long> &values) { return *std::min_element(values.cbegin(), values.cend()); }

/**
 * Arithmetic mean.
 * @param values
 * @return Arithmetic mean of the vector.
 */
inline long meanValue(const std::vector<long> &values) {
  return std::accumulate(values.cbegin(), values.cend(), 0l) / static_cast<long>(values.size());
}

/**
 * Median value.
 * @param values
 * @return Middle ((size-1) /2) of the sorted vector.
 */
template <class T>
inline T medianValue(std::vector<T> values) {
  if (values.empty()) return 0;

  std::sort(values.begin(), values.end());
  /// @todo C++20: replace by std::midpoint
  return values[(values.size() - 1) / 2];
}

/**
 * Optimal value according to passed strategy.
 * @param values
 * @param strategy For possible selector strategy choices see AutoPas::SelectorStrategy.
 * @return value or 0 if values is empty.
 */
inline long optimumValue(const std::vector<long> &values, SelectorStrategyOption strategy) {
  if (values.empty()) {
    return 0l;
  }

  switch (strategy) {
    case SelectorStrategyOption::fastestAbs: {
      return minValue(values);
    }
    case SelectorStrategyOption::fastestMean: {
      return meanValue(values);
    }
    case SelectorStrategyOption::fastestMedian: {
      return medianValue(values);
    }
      // no default to get appropriate compiler warnings
  }

  autopas::utils::ExceptionHandler::exception("OptionSelector: Unknown selector strategy {}!", strategy);

  return 0;
}

}  // namespace autopas::OptimumSelector
