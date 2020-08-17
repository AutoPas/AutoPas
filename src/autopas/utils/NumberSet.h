/**
 * @file NumberSet.h
 * @author Jan Nguyen
 * @date 09.06.19
 */

#pragma once

#include <set>
#include <sstream>
#include <vector>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Random.h"

namespace autopas {

/**
 * Virtual class describing a finite or infinite set of numbers
 */
template <class Number>
class NumberSet {
 public:
  virtual ~NumberSet() = default;

  /**
   * Create a copy of a NumberSet
   * @return
   */
  virtual std::unique_ptr<NumberSet> clone() const = 0;

  /**
   * Get a string representation of the set
   * @return string representation
   */
  virtual std::string to_string() const = 0;

  /**
   * Stream operator.
   * @param os
   * @param numberSet
   * @return
   */
  friend std::ostream &operator<<(std::ostream &os, const NumberSet &numberSet) {
    os << numberSet.to_string();
    return os;
  }

  /**
   * Indicates if the set is empty.
   * @return True if set is empty
   */
  virtual bool isEmpty() const = 0;

  /**
   * Indicates if the set is finite.
   * @return True if set is finite
   */
  virtual bool isFinite() const = 0;

  /**
   * Get size of set.
   * Only usable if set is finite.
   * @return size of set
   */
  virtual size_t size() const = 0;

  /**
   * Get the smallest number in the set
   * @return
   */
  virtual Number getMin() const = 0;
  /**
   * Get the largest number in the set
   * @return
   */
  virtual Number getMax() const = 0;

  /**
   * Get all numbers in the set. Only usable
   * if set is finite.
   * @return
   */
  virtual std::set<Number> getAll() const = 0;

  /**
   * Get a random number in the set.
   * @param rng random number generator
   * @return
   */
  virtual Number getRandom(Random &rng) const = 0;
  /**
   * Sample n points from the set. These points are
   * spaced evenly across the space.
   * @param n max samples
   * @param rng random number generator
   * @return samples randomly ordered vector of points
   */
  virtual std::vector<Number> uniformSample(size_t n, Random &rng) const = 0;

  /**
   * Sample up to n points from the set. These points are
   * spaced evenly across the space.
   * @param n max samples
   * @param rng random number generator
   * @return samples ordered set of points
   */
  virtual std::set<Number> uniformSampleSet(size_t n, Random &rng) const = 0;
};

}  // namespace autopas
