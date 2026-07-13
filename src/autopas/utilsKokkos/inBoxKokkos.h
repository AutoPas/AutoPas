/**
 * @file inBoxKokkos.h (Kokkos-clone of utils/inBox.h)
 * @date 09.07.2026
 * @author Luis Gall
 */

#pragma once

#include <Kokkos_Core.hpp>

namespace autopas::utilsKokkos {

// TODO: maybe find some way of merging inBox.h and inBoxKokkos.h via clever templating in the future
template <typename T>
KOKKOS_INLINE_FUNCTION
bool inBoxKokkos (const Kokkos::Array<T, 3>& position, const Kokkos::Array<T, 3>& low, const Kokkos::Array<T,3>& high) {
  static_assert(std::is_floating_point<T>::value, "inBoxKokkos assumes floating point types");

  bool inBox = true;
  for (int d = 0; d < 3; ++d) {
    const bool isLargerThanLower = position[d] >= low[d];
    const bool isSmallerThanHigher = position[d] < high[d];
    inBox = inBox and isLargerThanLower and isSmallerThanHigher;
  }
  return inBox;
}
}