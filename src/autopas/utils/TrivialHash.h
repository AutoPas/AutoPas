/**
 * @file TrivialHash.h
 * @author F. Gratl
 * @date 11/23/18
 */

#pragma once

// for std::size_t
#include <cstddef>

namespace autopas {

struct TrivialHash {
  template <typename T>
  std::size_t operator()(T t) const {
    return static_cast<std::size_t>(t);
  }
};

}  // namespace autopas
