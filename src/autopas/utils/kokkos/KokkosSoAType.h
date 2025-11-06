#pragma once
#include <tuple>

#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"
#include "autopas/utils/SoAType.h"

namespace autopas::utils::kokkos {

template <typename Space, typename... soatypes>
struct KokkosSoAType {};

template <typename Space, typename Particle_T, typename... soatypes>
struct KokkosSoAType<Space, autopas::utils::SoAType<Particle_T *, soatypes...>> {
  using Type = std::tuple<std::vector<Particle_T *>,
                          Kokkos::DualView<soatypes *, Kokkos::LayoutLeft, typename Space::memory_space>...>;
};

}  // namespace autopas::utils::kokkos
