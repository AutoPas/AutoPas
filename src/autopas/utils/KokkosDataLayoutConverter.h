/**
 *@file KokkosDataLayoutConverter.h
 *@author Luis Gall
 *@date 13.11.2025
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"

namespace autopas::utils {

template <class Particle_T>
class KokkosDataLayoutConverter {
 public:

  explicit KokkosDataLayoutConverter() {}

  template <class Input, class Output, std::size_t... I>
  void convertToSoA(Input &srcParticles, Output &dstParticles, size_t numParticles, std::index_sequence<I...> seq) {

    // AoS to SoA
    for (size_t i = 0; i < numParticles; ++i) {
      ((dstParticles. template operator()<I, false, true>(i) = srcParticles.getParticle(i).template get<static_cast<Particle_T::AttributeNames>(I+1)>()), ...); // I+1 as KokkosSoA does not contain ptr
    }

    dstParticles.template markAllModified<Kokkos::HostSpace::execution_space>(seq);
  }

  template <class Input, class Output, std::size_t... I>
  void convertToAoS(Input &srcParticles, Output &dstParticles, size_t numParticles, std::index_sequence<I...> seq) {

    srcParticles.template syncAll<Kokkos::HostSpace::execution_space>(seq);
    // SoA to AoS
    for (size_t i = 0; i < numParticles; ++i) {
      (dstParticles.getParticle(i). template set<static_cast<Particle_T::AttributeNames>(I+1)>(srcParticles. template operator()<I, false, true>(i)), ...); // I+1 as KokkosSoA does not contain ptr
    }
  }
};

}  // namespace autopas::utils
