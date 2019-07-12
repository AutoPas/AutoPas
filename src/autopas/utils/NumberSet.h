/**
 * @file NumberSet.h
 * @author Jan Nguyen
 * @date 09.06.19
 */

#pragma once

#include <set>
#include <sstream>
#include <vector>
#include "autopas/utils/ArrayUtils.h"
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
  virtual operator std::string() const = 0;

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
   * Sample n points from the set. These points are
   * spaced evenly across the space.
   * @param n max samples
   * @param rng random number generator
   * @return samples
   */
  virtual std::vector<Number> uniformSample(size_t n, Random &rng) const = 0;
};

/**
 * Class describing a finite set of numbers
 */
template <class Number>
class NumberSetFinite : public NumberSet<Number> {
 private:
  std::set<Number> _set;

 public:
  /**
   * Default Constructor: Empty set
   */
  NumberSetFinite() : _set() {}
  /**
   * Create a NumberSet of given values
   * @param values
   */
  NumberSetFinite(std::initializer_list<Number> values) : _set(values) {}
  /**
   * Create a NumberSet from a std::set
   * @param values
   */
  NumberSetFinite(std::set<Number> values) : _set(values) {}

  std::unique_ptr<NumberSet<Number>> clone() const override { return std::make_unique<NumberSetFinite>(*this); }

  operator std::string() const override { return "{" + ArrayUtils::to_string(_set) + "}"; }

  inline bool isEmpty() const override { return _set.empty(); }
  inline bool isFinite() const override { return true; }
  inline size_t size() const override { return _set.size(); }

  inline Number getMin() const override { return *_set.begin(); }
  inline Number getMax() const override { return *_set.rbegin(); }

  inline std::set<Number> getAll() const override { return _set; }

  std::vector<Number> uniformSample(size_t n, Random &rng) const override { return rng.uniformSample(_set, n); }
};

/**
 * Class describing an interval
 */
template <class Number>
class NumberInterval : public NumberSet<Number> {
 private:
  Number _min;
  Number _max;

 public:
  /**
   * Default Constructor: Create a range which only contains 0
   */
  NumberInterval() : _min(0.), _max(0.) {}
  /**
   * Create a range which only contains given value
   * @param val
   */
  NumberInterval(Number val) : _min(val), _max(val) {}
  /**
   * Create a range with given bounds
   * @param min Lower bound of range.
   * @param max Upper bound of range.
   */
  NumberInterval(Number min, Number max) : _min(min), _max(max) {
    if (min > max) {
      utils::ExceptionHandler::exception("NumberInterval: Invalid interval [{}, {}]", _min, _max);
    }
  }

  std::unique_ptr<NumberSet<Number>> clone() const override { return std::make_unique<NumberInterval>(*this); }

  operator std::string() const override {
    std::ostringstream ss;
    ss << "[" << _min << ", " << _max << "]";
    return ss.str();
  }

  inline bool isEmpty() const override { return false; }
  inline bool isFinite() const override { return _max == _min; }
  size_t size() const override {
    if (isFinite()) return 1ul;

    utils::ExceptionHandler::exception("NumberInterval.size: Interval is not finite [{}, {}]", _min, _max);
    return 0ul;
  }

  inline Number getMin() const override { return _min; }
  inline Number getMax() const override { return _max; }

  std::set<Number> getAll() const override {
    if (isFinite()) return {_min};

    utils::ExceptionHandler::exception("NumberInterval.getAll: Interval is not finite [{}, {}]", _min, _max);
    return {};
  }

  std::vector<Number> uniformSample(size_t n, Random &rng) const override {
    std::vector<Number> result;
    if (n == 0) {
      return result;
    }

    result.reserve(n);

    Number distance = (_max - _min) / (n - 1);
    for (size_t i = 0; i < (n - 1); ++i) {
      result.push_back(_min + distance * i);
    }
    // add max separatly, avoiding possible rounding errors
    result.push_back(_max);

    // randomize the sample
    rng.shuffle(std::begin(result), std::end(result));

    return result;
  }
};
}  // namespace autopas
