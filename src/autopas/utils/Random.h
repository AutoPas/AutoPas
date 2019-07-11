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

namespace autopas {

/**
 * Class for random algorithms.
 */
class Random {
 public:
  /**
   * Constructor
   */
  Random() : _rng(std::random_device()()) {}

  /**
   * Construct with seed.
   * @param seed
   */
  explicit Random(unsigned long seed) : _rng(seed) {}

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
   * @param pool
   * @param n number samples
   * @return samples
   */
  template <class T>
  std::vector<T> uniformSample(std::set<T> pool, size_t n) {
    std::vector<T> result;
    result.reserve(n + pool.size());

    // copy the whole set until result is full
    while (result.size() < n) {
      result.insert(result.end(), pool.begin(), pool.end());
    }

    // if too many elements added
    if (result.size() > n) {
      // randomize the last copy of the set
      size_t extra = result.size() - n;
      std::shuffle(std::end(result) - extra, std::end(result), _rng);

      // truncate the rest
      result.resize(n);
    }

    // randomize the sample
    std::shuffle(std::begin(result), std::end(result), _rng);

    return result;
  }

  /**
   * Reorders the elements in the given range [first, last) such that each possible
   * permutation of those elements has equal probability of appearance.
   * @param first
   * @param last
   */
  template <class RandomIt>
  inline void shuffle(RandomIt first, RandomIt last) {
    std::shuffle(first, last, _rng);
  }

 private:
  std::mt19937 _rng;
};

}  // namespace autopas
