/**
 * @file TrivialHash.h
 * @author F. Gratl
 * @date 11/23/18
 */

#pragma once

// for std::size_t
#include <cstddef>

namespace autopas {

/**
 * Trivial hash for enums.
 */
struct TrivialHash {
  /**
   * Trivial hash function
   * @tparam T
   * @param t
   * @return Integer representation of t
   */
  template <typename T>
  std::size_t operator()(T t) const {
    return static_cast<std::size_t>(t);
  }
};

}  // namespace autopas
