#pragma once

#include "Kokkos_Core.hpp"
#include "autopas/utils/SoAStorage.h"

namespace autopas::utils::kokkos {

template <class SoAArraysType>
class KokkosSoA {
 public:
  KokkosSoA() = default;
  KokkosSoA(const KokkosSoA &soa) = default;

  void resize(size_t newSize) {
    _soaStorage.apply([newSize](auto view) { Kokkos::resize(view, newSize); });
  }
  [[nodiscard]] constexpr size_t size() const { return _soaStorage.template get<0>().extent(0); }

 private:
  SoAStorage<SoAArraysType> _soaStorage;
};

}  // namespace autopas::utils::kokkos