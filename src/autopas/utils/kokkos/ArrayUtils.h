#pragma once

#include "Kokkos_Core.hpp"

namespace autopas::utils::kokkos::ArrayUtils {
template <class ScalarType, int N>
struct array_type {
  std::array<ScalarType, N> the_array;

  KOKKOS_INLINE_FUNCTION
  array_type() {
    for (int i = 0; i < N; i++) {
      the_array[i] = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  array_type(const array_type &rhs) {
    for (int i = 0; i < N; i++) {
      the_array[i] = rhs.the_array[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  ScalarType &operator[](int i) { return the_array[i]; }

  KOKKOS_INLINE_FUNCTION
  const ScalarType &operator[](int i) const { return the_array[i]; }

  KOKKOS_INLINE_FUNCTION
  array_type &operator+=(const array_type &src) {
    for (int i = 0; i < N; i++) {
      the_array[i] += src.the_array[i];
    }
    return *this;
  }
};

typedef array_type<double, 3> Vector3;
}  // namespace array_utils

namespace Kokkos {
template <>
struct reduction_identity<autopas::utils::kokkos::ArrayUtils::Vector3> {
  KOKKOS_FORCEINLINE_FUNCTION static autopas::utils::kokkos::ArrayUtils::Vector3 sum() { return {}; }
};

}  // namespace Kokkos
