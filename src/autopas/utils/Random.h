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
  explicit Random(unsigned long seed = std::random_device()()) : std::mt19937(seed) {}

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
   * Get a uniformly random object from the given set.
   * @param pool set
   * @return random element
   */
  template <class T>
  T pickRandom(std::set<T> pool) {
    std::uniform_int_distribution<size_t> distr(0ul, pool.size() - 1ul);
    size_t pos = distr(*this);

    auto it = pool.begin();
    for (size_t i = 0; i < pos; ++i) ++it;

    return *it;
  }
};

}  // namespace autopas
