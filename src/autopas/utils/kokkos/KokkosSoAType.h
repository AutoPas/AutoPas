#pragma once
#include <tuple>

#include "KokkosSpace.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"
#include "autopas/utils/SoAType.h"

namespace autopas::utils::kokkos {

template <typename... soatypes>
struct KokkosSoAType {};

template <typename Particle_T, typename... soatypes>
struct KokkosSoAType<autopas::utils::SoAType<Particle_T *, soatypes...>> {
  using Type = std::tuple<Kokkos::DualView<uintptr_t * , Kokkos::LayoutLeft, DeviceSpace>, Kokkos::DualView<soatypes *, Kokkos::LayoutLeft, DeviceSpace>...>;
};

}  // namespace autopas::utils::kokkos
