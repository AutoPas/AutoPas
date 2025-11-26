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

  explicit KokkosDataLayoutConverter(DataLayoutOption dataLayout)
      : _dataLayout(dataLayout) {}

  template <class Input, class Output, std::size_t... I>
  void loadDataLayout(Input &srcParticles, Output &dstParticles, size_t numParticles, std::index_sequence<I...>) {

    // AoS to SoA
    for (size_t i = 0; i < numParticles; ++i) {
      (dstParticles. template set<I>(srcParticles(i).template get<static_cast<Particle_T::AttributeNames>(I+1)>(), i), ...); // I+1 as KokkosSoA does not contain ptr
    }
  }

  template <class Input, class Output, std::size_t... I>
  void storeDataLayout(Input &srcParticles, Output &dstParticles, size_t numParticles, std::index_sequence<I...>) {

    // SoA to AoS
    for (size_t i = 0; i < numParticles; ++i) {
      (dstParticles(i). template set<static_cast<Particle_T::AttributeNames>(I+1)>(srcParticles. template get<I>(i)), ...); // I+1 as KokkosSoA does not contain ptr
    }
  }

 private:

  DataLayoutOption _dataLayout;
};

}  // namespace autopas::utils
