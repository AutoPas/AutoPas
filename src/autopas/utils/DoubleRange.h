/**
 * @file DoubleRange.h
 * @author Jan Nguyen
 * @date 17.05.19
 */

#pragma once

#include "autopas/autopasIncludes.h"

namespace autopas {

/**
 * Class describing a range of doubles
 */
class DoubleRange {
 private:
  double _min;
  double _max;

 public:
  /**
   * Default Constructor: Create a range which only contains 0
   */
  DoubleRange() : DoubleRange(0., 0.) {}
  /**
   * Create a range which only contains given value
   * @param val
   */
  DoubleRange(double val) : _min(val), _max(val) {}
  /**
   * Create a range with given bounds
   * @param min Lower bound of range.
   * @param max Upper bound of range.
   */
  DoubleRange(double min, double max) : _min(min), _max(max) {}

  /**
   * Gets lower bound of range
   * @return lower bound
   */
  double getMin() const { return _min; }
  /**
   * Gets upper bound of range
   * @return upper bound
   */
  double getMax() const { return _max; }

  /**
   * Create n equidistant values with getMin()
   * being the smallest and getMax() the greatest.
   * @param n
   * @return
   */
  std::vector<double> range(unsigned n) const {
    std::vector<double> result;
    result.reserve(n);

    double distance = (_max - _min) / (n - 1);
    for (unsigned i = 0; i < n; ++i) {
      result.push_back(_min + distance * i);
    }

    return result;
  }
};
}  // namespace autopas
