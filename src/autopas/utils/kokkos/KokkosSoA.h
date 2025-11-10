#pragma once

#include "KokkosSoAType.h"
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
  [[nodiscard]] constexpr size_t size() const { return _soaStorage.template get<0>().size(); }

  void clear() {
    _soaStorage.apply([](auto &view) { Kokkos::resize(view, 0); });
  }
  template <size_t AttributeName, typename Space>
  inline constexpr void markModified() {
    auto &view = _soaStorage.template get<AttributeName>();
    view.template modify<Space>();
  }
  template <size_t AttributeName, typename Space>
  inline constexpr void sync() {
    auto &view = _soaStorage.template get<AttributeName>();
    view.template sync<Space>();
  }

  template <typename Space>
  inline constexpr void syncAll() {
    std::apply([](auto &view) { view.template sync<Space>(); });
  }

  template <typename Space>
  inline constexpr void markModifiedAll() {
    std::apply([](auto &view) { view.template modify<Space>(); });
  }

 private:
  SoAStorage<typename KokkosSoAType<SoAArraysType>::Type> _soaStorage;
};

}  // namespace autopas::utils::kokkos