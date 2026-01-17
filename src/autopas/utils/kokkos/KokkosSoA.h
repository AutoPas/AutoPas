#pragma once

#include "KokkosSoAType.h"
#include "Kokkos_Core.hpp"
#include "autopas/utils/SoAStorage.h"

namespace autopas::utils::kokkos {

template <class SoAArraysType>
class KokkosSoA {
 public:
  KokkosSoA() {
    _soaStorage.apply([](auto &view) { view = typename std::decay_t<decltype(view)>("Label", 0); });
  };
  KokkosSoA(const KokkosSoA &soa) = default;

  void resize(size_t newSize) {
    Kokkos::fence();
    _soaStorage.apply([newSize](auto &view) { view.resize(newSize); });
    syncAll<DeviceSpace>();
    Kokkos::fence();
    syncAll<HostSpace>();
    Kokkos::fence();
  }
  [[nodiscard]] constexpr size_t size() const { return _soaStorage.template get<0>().extent(0); }

  void clear() {
    _soaStorage.apply([](auto &view) { view.resize(0); });
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
    _soaStorage.apply([](auto &view) {
      view.template sync<Space>();
    });
  }

  template <typename Space>
  inline constexpr void markModifiedAll() {
    _soaStorage.apply([](auto &view) {
      view.template modify<Space>();
    });
  }

  template <size_t AttributeName>
  inline constexpr auto &get() const {
    return _soaStorage.template get<AttributeName>();
  }

  template <typename Space, typename PermView>
  inline void sort(Space space, PermView permutation) {
    _soaStorage.apply([&](auto &view) {
      apply_permutation(space, permutation, view.d_view);
      view.template modify<Space>();
    });
  }

  template <typename Space, typename View, typename PermView>
  void apply_permutation(Space space, PermView permutations, View &view) {
    auto temp = Kokkos::create_mirror(Kokkos::WithoutInitializing, space, view);
    auto policy = Kokkos::RangePolicy<Space>(0, size());

    Kokkos::parallel_for(
        "apply permutation", policy, KOKKOS_LAMBDA(const int index) {
          // This is a hack to get it working. I have to declare them, because otherwise I can't use constexpr
          // without constexpr I can't handle the dimensionality
          // And they are not allowed to be shadowed
          // View are copied shallow
          View temp1 = temp;
          View view1 = view;
          auto new_index = permutations(index);
          if constexpr (View::rank == 1) {
            // uncoalesced memory access, but unavoidable at this point
            temp1(index) = view1(new_index);
          } else if constexpr (View::rank == 2) {
            temp1(index, 0) = view1(new_index, 0);
            temp1(index, 1) = view1(new_index, 1);
            temp1(index, 2) = view1(new_index, 2);
          }
        });
    view = temp;
  };

 private:
  SoAStorage<typename KokkosSoAType<SoAArraysType>::Type> _soaStorage;
};

}  // namespace autopas::utils::kokkos