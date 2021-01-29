/**
 * @file Random.h
 *
 * @date 19.06.19
 * @author Jan Nguyen
 */

#pragma once

#include <algorithm>
#include <random>
#include <set>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Class for random algorithms.
 */
class Random : public std::mt19937 {
 public:
  /**
   * Constructor
   * @param seed
   */
  explicit Random(uint64_t seed = std::random_device()()) : std::mt19937(seed) {}

  /**
   * Class should not be copied constructed
   * @param other
   */
  Random(const Random &other) = delete;

  /**
   * Class should not be copied assigned
   * @param other
   * @return
   */
  Random &operator=(const Random &other) = delete;

  /**
   * Sample n points from the pool. Each element in the pool will
   * appear about the same number of times in the sample.
   * @tparam T type of the elements
   * @param poolBegin
   * @param poolEnd
   * @param n number samples
   * @return samples
   */
  template <class Iter>
  std::vector<typename std::iterator_traits<Iter>::value_type> uniformSample(Iter poolBegin, Iter poolEnd, size_t n) {
    if (poolBegin == poolEnd) {
      autopas::utils::ExceptionHandler::exception("Random.uniformSample: Cannot sample from empty set.");
    }

    std::vector<typename std::iterator_traits<Iter>::value_type> result;
    result.reserve(n);

    // copy the whole set until result is full
    while (result.size() < n) {
      result.insert(std::end(result), poolBegin, poolEnd);
    }

    // if too many elements added
    if (result.size() > n) {
      // randomize the last copy of the set
      size_t extra = result.size() - n;
      std::shuffle(std::end(result) - extra, std::end(result), *this);

      // truncate the rest
      result.resize(n);
    }

    // randomize the sample
    std::shuffle(std::begin(result), std::end(result), *this);

    return result;
  }

  /**
   * Sample n points from the set {min;min+1;...;max}. Each element in the set will
   * appear about the same number of times in the sample.
   * @param min
   * @param max
   * @param n number samples
   * @return samples
   */
  std::vector<size_t> uniformSample(size_t min, size_t max, size_t n) {
    std::set<size_t> allAllowed;
    for (size_t i = min; i <= max; ++i) {
      allAllowed.insert(i);
    }
    return uniformSample(allAllowed.begin(), allAllowed.end(), n);
  }

  /**
   * Get a uniformly random object from the given set.
   * @param pool set
   * @return random element
   */
  template <class T>
  T pickRandom(std::set<T> pool) {
    std::uniform_int_distribution<size_t> distr(0ul, pool.size() - 1ul);
    size_t pos = distr(*this);

    auto it = pool.begin();
    std::advance(it, pos);

    return *it;
  }

  /**
   * Pick up to n random elements from the set.
   * Each elements has the same probability to be chosen
   * @param pool set
   * @param n number of elements
   * @return n random elements
   */
  template <class T>
  std::set<T> randomSubset(std::set<T> pool, size_t n) {
    size_t size = std::min(n, pool.size());

    // create randomly shuffled vector of points to all elements in the pool
    std::vector<const T *> pointerVec;
    pointerVec.reserve(pool.size());
    for (const T &element : pool) {
      pointerVec.push_back(&element);
    }
    std::shuffle(pointerVec.begin(), pointerVec.end(), *this);

    // return first elements of the shuffled vector
    std::set<T> result;
    for (size_t i = 0; i < size; ++i) {
      result.insert(*pointerVec[i]);
    }

    return result;
  }
};

}  // namespace autopas
