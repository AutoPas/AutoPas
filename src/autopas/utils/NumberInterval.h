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
    // add max separatly, avoiding possible rounding errors
    result.push_back(_max);

    // randomize the sample
    std::shuffle(std::begin(result), std::end(result), rng);

    return result;
  }
};
}  // namespace autopas
