/**
 *@file KokkosDataLayoutConverter.h
 *@author Luis Gall
 *@date 13.11.2025
 */

#pragma once

namespace autopas::utilsKokkos {

class KokkosDataLayoutConverter {
 public:

  template <typename ExecSpace, bool useHostView, class Input, class Output, std::size_t... I>
  static void convertToSoA(Input &srcParticles, Output &dstParticles, size_t numParticles, std::index_sequence<I...> seq) {
    // AoS to SoA
    Kokkos::parallel_for("autopas::KokkosDataLayoutConverter::convertToSoA", Kokkos::RangePolicy<ExecSpace>(0, numParticles), KOKKOS_LAMBDA (size_t i) {
      ((dstParticles. template operator()<I, useHostView>(i) = srcParticles.template operator()<I+1, useHostView>(i)), ...); // I+1 as KokkosSoA does not contain ptr
    });
  }

  template <typename ExecSpace, bool useHostView, class Input, class Output, std::size_t... I>
  static void convertToAoS(Input &srcParticles, Output &dstParticles, size_t numParticles, std::index_sequence<I...> seq) {
    // SoA to AoS
    Kokkos::parallel_for("autopas::KokkosDataLayoutConverter::convertToAoS", Kokkos::RangePolicy<ExecSpace>(0,numParticles), KOKKOS_LAMBDA (size_t i) {
      ((dstParticles.template operator()<I+1, useHostView>(i) = srcParticles.template operator()<I, useHostView>(i)), ...); // I+1 as KokkosSoA does not contain ptr
    });
  }
};

}  // namespace autopas::utils
