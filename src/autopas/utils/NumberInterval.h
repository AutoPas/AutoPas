/**
 * @file NumberInterval.h
 * @author F. Gratl
 * @date 11/15/19
 */

#pragma once
#include "NumberSet.h"

namespace autopas {
/**
 * Class describing an interval
 */
template <class Number>
class NumberInterval : public NumberSet<Number> {
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

  /**
   * Setter for NumberInterval
   * @param numbers One or two values, like the available constructors for NumberInterval.
   * If two are provided the smaller one is assumed to be the min value.
   */
  inline void resetValues(std::set<Number> &numbers) override {
    if (numbers.size() == 1) {
      _min = *numbers.begin();
      _max = *numbers.begin();
    } else if (numbers.size() == 2) {
      _min = std::min(*numbers.begin(), *++numbers.begin());
      _max = std::max(*numbers.begin(), *++numbers.begin());

    } else {
      utils::ExceptionHandler::exception(
          "NumberInterval::resetValues: NumberInterval constructor takes exactly"
          " one or two values");
    }
  }

  std::string to_string() const override {
    std::ostringstream ss;
    ss << _min << "-" << _max;
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

  inline Number getRandom(Random &rng) const override {
    std::uniform_real_distribution<Number> distr(_min, _max);
    return distr(rng);
  }
  std::vector<Number> uniformSample(size_t n, Random &rng) const override {
    std::vector<Number> result;
    if (n == 0) {
      return result;
    } else if (n == 1) {
      // if only one sample choose middle
      result.reserve(1);
      result.push_back((_max + _min) / 2);
      return result;
    }

    result.reserve(n);

    Number distance = (_max - _min) / (n - 1);
    for (size_t i = 0; i < (n - 1); ++i) {
      result.push_back(_min + distance * i);
    }
    // add max separately, avoiding possible rounding errors
    result.push_back(_max);

    // randomize the sample
    std::shuffle(std::begin(result), std::end(result), rng);

    return result;
  }

  std::set<Number> uniformSampleSet(size_t n, Random &rng) const override {
    std::set<Number> result;
    if (n == 0) {
      return result;
    } else if (isFinite()) {
      result.insert(_min);
      return result;
    } else if (n == 1) {
      // if only one sample choose median
      result.insert(getMedian());
      return result;
    }

    Number distance = (_max - _min) / (n - 1);
    for (size_t i = 0; i < (n - 1); ++i) {
      result.insert(_min + distance * i);
    }
    // add max separately, avoiding possible rounding errors
    result.insert(_max);

    return result;
  }

  Number getMedian() const override { return (_max + _min) / 2; }

 private:
  Number _min;
  Number _max;
};
}  // namespace autopas
