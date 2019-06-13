/**
 * @file NumberSet.h
 * @author Jan Nguyen
 * @date 09.06.19
 */

#pragma once

#include <random>
#include <set>
#include <sstream>
#include <vector>
#include "autopas/utils/ExceptionHandler.h"

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
  virtual std::vector<Number> uniformSample(unsigned n, std::default_random_engine& rng) const = 0;
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

  operator std::string() const override {
    auto it = _set.begin();
    std::stringstream ss;
    ss << "{" << *it;
    for (++it; it != _set.end(); ++it) {
      ss << ", " << *it;
    }
    ss << "}";

    return ss.str();
  }

  bool isEmpty() const override { return _set.empty(); }
  bool isFinite() const override { return true; }

  Number getMin() const override { return *_set.begin(); }
  Number getMax() const override { return *_set.rbegin(); }

  std::set<Number> getAll() const override { return _set; }

  std::vector<Number> uniformSample(unsigned n, std::default_random_engine& rng) const override {
    std::vector<Number> result;
    result.reserve(n + _set.size());

    // copy the whole set until result is full
    while (result.size() < n) {
      for (auto val : _set) {
        result.push_back(val);
      }
    }

    // if too many elements added
    if (result.size() > n) {
      // randomize the last copy of the set
      size_t extra = result.size() - n;
      std::shuffle(std::end(result) - extra, std::end(result), rng);

      // truncate the rest
      result.resize(n);
    }

    // randomize the sample
    std::shuffle(std::begin(result), std::end(result), rng);

    return result;
  }
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
    std::stringstream ss;
    ss << "[" << _min << ", " << _max << "]";
    return ss.str();
  }

  Number getMin() const override { return _min; }
  Number getMax() const override { return _max; }

  bool isEmpty() const override { return false; }
  bool isFinite() const override { return _max == _min; }

  std::set<Number> getAll() const override {
    if (isFinite()) return {_min};

    utils::ExceptionHandler::exception("NumberInterval: Interval is not finite [{}, {}]", _min, _max);
    return {};
  }

  std::vector<Number> uniformSample(unsigned n, std::default_random_engine& rng) const override {
    std::vector<Number> result;
    if (n == 0) {
      return result;
    }

    result.reserve(n);

    Number distance = (_max - _min) / (n - 1);
    for (unsigned i = 0; i < (n - 1); ++i) {
      result.push_back(_min + distance * i);
    }
    // add max separatly, avoiding possible rounding errors
    result.push_back(_max);

    // randomize the sample
    std::shuffle(std::begin(result), std::end(result), rng);

    return result;
  }
};
}  // namespace autopas
