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

  explicit DSKokkosTraversalInterface(DataLayoutOption dataLayout) : _dataLayout(dataLayout), _converter(_dataLayout) {};

  virtual ~DSKokkosTraversalInterface() = default;

  void setOwnedAoSToTraverse(utils::KokkosAoS<HostSpace, Particle_T>& inputAoS) {
    // Target Layout is SoA
    const size_t numParticles = inputAoS.size();

    // TODO: this should in the end become a mirror view for avoiding data copies where they are not necessary (i.e. HostSpace = DeviceSpace)
    KokkosSoAType<HostSpace> hostOwnedParticlesSoA {numParticles};
    _deviceOwnedParticlesSoA.resize(numParticles);

    constexpr size_t tupleSize = _deviceOwnedParticlesSoA.tupleSize();
    constexpr auto I = std::make_index_sequence<tupleSize>();

    if (_dataLayout == DataLayoutOption::soa) {
      _converter.convertToSoA(inputAoS, hostOwnedParticlesSoA, numParticles, I);
    }
    else if (_dataLayout == DataLayoutOption::aos) {
      // No Op as no conversion is required
    }

    _deviceOwnedParticlesSoA.copyFrom(hostOwnedParticlesSoA, I);
  }

  void setOwnedSoAToTraverse(KokkosSoAType<HostSpace>& inputSoA) {

    const size_t numParticles = inputSoA.size();

    utils::KokkosAoS<HostSpace, Particle_T> hostOwnedParticlesAoS {numParticles};
    Kokkos::realloc(_ownedParticlesAoS, numParticles);

    // Target Layout is SoA
    if (_dataLayout == DataLayoutOption::soa) {
      // No Op as no conversion is required
    }
    else if (_dataLayout == DataLayoutOption::aos) {
      constexpr size_t tupleSize = inputSoA.tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      _converter.convertToAoS(inputSoA, hostOwnedParticlesAoS, numParticles, I);
    }

    Kokkos::deep_copy(_ownedParticlesAoS, hostOwnedParticlesAoS.getView());
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

private:

  DataLayoutOption _dataLayout;
  utils::KokkosDataLayoutConverter<Particle_T> _converter;
};
}