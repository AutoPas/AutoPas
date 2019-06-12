/**
 * @file DoubleSet.h
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
 * Virtual class describing a finite or infinite set of doubles
 */
class DoubleSet {
 public:
  virtual ~DoubleSet() = default;

  /**
   * Create a copy of a DoubleSet
   * @return
   */
  virtual std::unique_ptr<DoubleSet> clone() const = 0;

  /**
   * Get string representation of set
   */
  virtual operator std::string() const = 0;

  /**
   * Checks if the set is finite.
   * @return True if set is finite
   */
  virtual bool isFinite() const = 0;

  /**
   * Get the smallest double in the set
   * @return
   */
  virtual double getMin() const = 0;
  /**
   * Get the largest double in the set
   * @return
   */
  virtual double getMax() const = 0;

  /**
   * Get number of doubles in the set.
   * Only usable if set is finite.
   * @return size of set
   */
  virtual size_t size() const = 0;

  /**
   * Get all doubles in the set. Only usable
   * if set is finite.
   * @return set as std::set
   */
  virtual std::set<double> getAll() const = 0;

  /**
   * Sample n points from the set. These points are
   * spaced evenly across the space.
   * @param n max samples
   * @param rng random number generator
   * @return samples
   */
  virtual std::vector<double> uniformSample(unsigned n, std::default_random_engine& rng) const = 0;
};

/**
 * Class describing a finite set of doubles
 */
class DoubleFiniteSet : public DoubleSet {
 private:
  std::set<double> _set;

 public:
  /**
   * Default Constructor: Empty set
   */
  DoubleFiniteSet() : _set() {}
  /**
   * Create a DoubleSet of given values
   * @param values
   */
  DoubleFiniteSet(std::initializer_list<double> values) : _set(values) {}
  /**
   * Create a DoubleSet from a std::set
   * @param values
   */
  DoubleFiniteSet(std::set<double> values) : _set(values) {}

  std::unique_ptr<DoubleSet> clone() const override { return std::make_unique<DoubleFiniteSet>(*this); }

  operator std::string() const override {
    if (_set.empty()) {
      return "{ }";
    }

    auto it = _set.begin();
    std::stringstream ss;
    ss << "{" << *it;
    for (++it; it != _set.end(); ++it) {
      ss << ", " << *it;
    }
    ss << "}";

    return ss.str();
  }

  bool isFinite() const override { return true; }

  double getMin() const override { return *_set.begin(); }
  double getMax() const override { return *_set.rbegin(); }

  size_t size() const override { return _set.size(); }

  std::set<double> getAll() const override { return _set; }

  std::vector<double> uniformSample(unsigned n, std::default_random_engine& rng) const override {
    std::vector<double> result;
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
 * Class describing a interval
 */
class DoubleInterval : public DoubleSet {
 private:
  double _min;
  double _max;

 public:
  /**
   * Default Constructor: Create a range which only contains 0
   */
  DoubleInterval() : _min(0.), _max(0.) {}
  /**
   * Create a range which only contains given value
   * @param val
   */
  DoubleInterval(double val) : _min(val), _max(val) {}
  /**
   * Create a range with given bounds
   * @param min Lower bound of range.
   * @param max Upper bound of range.
   */
  DoubleInterval(double min, double max) : _min(min), _max(max) { assert(min <= max); }

  std::unique_ptr<DoubleSet> clone() const override { return std::make_unique<DoubleInterval>(*this); }

  operator std::string() const override {
    std::stringstream ss;
    ss << "[" << _min << ", " << _max << "]";
    return ss.str();
  }

  bool isFinite() const override { return _max == _min; }

  double getMin() const override { return _min; }
  double getMax() const override { return _max; }

  size_t size() const override {
    if (isFinite()) return 1;

    utils::ExceptionHandler::exception("DoubleInterval.getSize: Interval is not finite [{}, {}]", _min, _max);
    return 0;
  }

  std::set<double> getAll() const override {
    if (isFinite()) return {_min};

    utils::ExceptionHandler::exception("DoubleInterval.getAll: Interval is not finite [{}, {}]", _min, _max);
    return {};
  }

  std::vector<double> uniformSample(unsigned n, std::default_random_engine& rng) const override {
    std::vector<double> result;
    result.reserve(n);

    double distance = (_max - _min) / (n - 1);
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
