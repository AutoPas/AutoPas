/**
* @file DSKokkosTraversalInterface.h
 * @date 05.11.2025
 * @author Luis Gall
 */

#pragma once

#include <Kokkos_Core.hpp>
#include "autopas/utils/KokkosAoS.h"
#include "autopas/utils/KokkosDataLayoutConverter.h"
#include "autopas/options/DataLayoutOption.h"

#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
using DeviceSpace = Kokkos::HostSpace;
#endif

using HostSpace = Kokkos::HostSpace;

namespace autopas {
template <class Particle_T>
class DSKokkosTraversalInterface {

public:

  template <typename Space>
  using KokkosSoAType = Particle_T::template KokkosSoAArraysType<Space>;

  explicit DSKokkosTraversalInterface() {};

  virtual ~DSKokkosTraversalInterface() = default;

  void setOwnedAoSToTraverse(utils::KokkosAoS<HostSpace, Particle_T>& inputAoS) {
    Kokkos::realloc(_ownedParticlesAoS, inputAoS.size());
    Kokkos::deep_copy(_ownedParticlesAoS, inputAoS.getView());
  }

  void setOwnedSoAToTraverse(KokkosSoAType<HostSpace>& inputSoA) {
    constexpr size_t tupleSize = _deviceOwnedParticlesSoA.tupleSize();
    constexpr auto I = std::make_index_sequence<tupleSize>();

    _deviceOwnedParticlesSoA.resize(inputSoA.size());
    _deviceOwnedParticlesSoA.copyFrom(inputSoA, I);
  }

  void retrieveOwnedAoS(utils::KokkosAoS<HostSpace, Particle_T>& outputAoS) {
    Kokkos::deep_copy(outputAoS.getView(), _ownedParticlesAoS);
  }

  void retrieveOwnedSoA(KokkosSoAType<HostSpace>& outputSoA) {
    constexpr size_t tupleSize = _deviceOwnedParticlesSoA.tupleSize();
    constexpr auto I = std::make_index_sequence<tupleSize>();

    outputSoA.copyFrom(_deviceOwnedParticlesSoA, I);
  }

  template <class IncomingParticles_T>
  void setHaloToTraverse(IncomingParticles_T &particles) {
    // TODO
  };

protected:

  Kokkos::View<Particle_T*, DeviceSpace> _ownedParticlesAoS {0};
  Kokkos::View<Particle_T*, DeviceSpace> _haloParticlesAoS {0};

  KokkosSoAType<DeviceSpace> _deviceOwnedParticlesSoA {0};
  KokkosSoAType<DeviceSpace> _deviceHaloParticlesSoA {0};
};
}