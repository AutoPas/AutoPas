/**
 * @file Random.h
 *
 * @date 19.06.19
 * @author Jan Nguyen
 */

#pragma once

#include <algorithm>
#include <iterator>
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
  std::vector<typename std::iterator_traits<Iter>::value_type> uniformSample(Iter poolBegin, Iter poolEnd, std::size_t n) {
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
      std::size_t extra = result.size() - n;
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
  std::vector<std::size_t> uniformSample(std::size_t min, std::size_t max, std::size_t n) {
    std::set<std::size_t> allAllowed;
    for (std::size_t i = min; i <= max; ++i) {
      allAllowed.insert(i);
    }
    return uniformSample(allAllowed.begin(), allAllowed.end(), n);
  }

  /**
   * Get a uniformly randomly selected object from the given container.
   *
   * @tparam Container Type of the container. Must support std::begin().
   * @tparam Elem Type of the elements in the container.
   * @param pool Container from which to select an element.
   * @return Randomly selected element.
   */
  template <class Container>
  auto pickRandom(const Container &pool) {
    std::uniform_int_distribution<std::size_t> distr(0ul, pool.size() - 1ul);
    std::size_t pos = distr(*this);

    auto it = std::begin(pool);
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
  std::set<T> randomSubset(std::set<T> pool, std::size_t n) {
    std::size_t size = std::min(n, pool.size());

    // create randomly shuffled vector of points to all elements in the pool
    std::vector<const T *> pointerVec;
    pointerVec.reserve(pool.size());
    for (const T &element : pool) {
      pointerVec.push_back(&element);
    }
    std::shuffle(pointerVec.begin(), pointerVec.end(), *this);

    // return first elements of the shuffled vector
    std::set<T> result;
    for (std::size_t i = 0; i < size; ++i) {
      result.insert(*pointerVec[i]);
    }

    return result;
  }
};

}  // namespace autopas
