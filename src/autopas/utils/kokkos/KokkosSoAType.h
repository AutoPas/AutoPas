#pragma once
#include <tuple>

#include "Kokkos_Core.hpp"
#include "autopas/utils/SoAType.h"

namespace autopas::utils::kokkos {

template <typename Space, typename Layout, typename... soatypes>
struct KokkosSoAType {
  using Type = std::tuple<Kokkos::View<soatypes *, Layout, typename Space::memory_space>...>;
};

//
// template <typename Space, typename Layout, typename... soatypes>
// struct KokkosSoAType<Space, Layout, autopas::utils::SoAType<soatypes...>> {
//   using Type = std::tuple<Kokkos::View<soatypes *, Layout, typename Space::memory_space>...>;
// };

// // Partial template specialization on the SoAType::Type,
// template <typename Space, typename Layout, typename... soatypes>
// struct KokkosSoAType<Space, Layout, std::tuple<std::vector<soatypes, AlignedAllocator<soatypes>>...>> {
//   using Type = std::tuple<Kokkos::View<soatypes *, Layout, typename Space::memory_space>...>;
// };

// Partial template specialization on the SoAType::Type, Pointer to the first parameter is skipped, because it's a
// pointer to the particle
template <typename Space, typename Layout, typename Particle_T, typename... soatypes>
struct KokkosSoAType<Space, Layout,
                     std::tuple<std::vector<Particle_T *, AlignedAllocator<Particle_T *>>,
                                std::vector<soatypes, AlignedAllocator<soatypes>>...>> {
  using Type = std::tuple<Kokkos::View<soatypes *, Layout, typename Space::memory_space>...>;
};

}  // namespace autopas::utils::kokkos
