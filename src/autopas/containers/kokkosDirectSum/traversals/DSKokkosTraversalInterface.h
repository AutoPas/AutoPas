/**
* @file DSKokkosTraversalInterface.h
 * @date 05.11.2025
 * @author Luis Gall
 */

#pragma once

#include <Kokkos_Core.hpp>
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

  explicit DSKokkosTraversalInterface(DataLayoutOption dataLayout) : _dataLayout(dataLayout), _converter(_dataLayout) {};

  virtual ~DSKokkosTraversalInterface() = default;

  virtual void setOwnedToTraverse(Kokkos::View<Particle_T*, HostSpace> &particles) {

    // TODO: consider chosen data layout
    const size_t numParticles = particles.extent(0);
    Kokkos::realloc(_ownedParticlesAoS, numParticles);

    Kokkos::deep_copy(_ownedParticlesAoS, particles);

    KokkosSoAType<HostSpace> hostOwnedParticlesSoA {};
    _deviceOwnedParticleSoA.resize(numParticles);
    hostOwnedParticlesSoA.resize(numParticles);

    constexpr size_t tupleSize = hostOwnedParticlesSoA.tupleSize();
    constexpr auto I = std::make_index_sequence<tupleSize>();

    _converter.loadDataLayout(particles, hostOwnedParticlesSoA, numParticles, I);
    _deviceOwnedParticleSoA.copyFrom(hostOwnedParticlesSoA, I);
  }

  virtual void setHaloToTraverse(const Kokkos::View<Particle_T*, HostSpace> &particles) {
    // TODO
  };

protected:

  Kokkos::View<Particle_T*, DeviceSpace> _ownedParticlesAoS;
  Kokkos::View<Particle_T*, DeviceSpace> _haloParticlesAoS;

  template <typename Space>
  using KokkosSoAType = Particle_T::template KokkosSoAArraysType<Space>;

  KokkosSoAType<DeviceSpace> _deviceOwnedParticleSoA;

private:
  DataLayoutOption _dataLayout;
  utils::KokkosDataLayoutConverter<Particle_T> _converter;
};
}