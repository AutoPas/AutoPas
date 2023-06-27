/**
 * @file NumberSetFinite.h
 * @author F. Gratl
 * @date 11/15/19
 */

#pragma once

#include <set>

#include "ArrayUtils.h"
#include "NumberSet.h"

namespace autopas {
/**
 * Class describing a finite set of numbers
 */
template <class Number>
class NumberSetFinite final : public NumberSet<Number> {
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
  /**
   * Setter for NumberSetFinite
   * @param numbers The set of numbers the new NumberSetFinite represents
   */
  inline void resetValues(const std::set<Number> &numbers) override { _set = numbers; }

  std::string to_string() const override { return "" + utils::ArrayUtils::to_string(_set) + ""; }

  inline bool isEmpty() const override { return _set.empty(); }
  inline bool isFinite() const override { return true; }
  inline size_t size() const override { return _set.size(); }

  inline Number getMin() const override { return *_set.begin(); }
  inline Number getMax() const override { return *_set.rbegin(); }

  inline std::set<Number> getAll() const override { return _set; }

  inline Number getRandom(Random &rng) const override { return rng.pickRandom(_set); }
  std::vector<Number> uniformSample(size_t n, Random &rng) const override {
    return rng.uniformSample(_set.begin(), _set.end(), n);
  }
  std::set<Number> uniformSampleSet(size_t n, Random &rng) const override { return rng.randomSubset(_set, n); }
  Number getMedian() const override {
    auto it = std::begin(_set);
    std::advance(it, _set.size() / 2);
    return *it;
  }

 private:
  /**
   * Set containing all numbers.
   */
  std::set<Number> _set;
};
}  // namespace autopas
