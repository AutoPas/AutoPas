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

    KokkosSoAType<HostSpace> ownedSoA {particles.extent(0)};


    // TODO: decide where data conversion should take place, GPU/CPU
    auto ownedView = Kokkos::create_mirror_view(DeviceSpace{}, particles);
    Kokkos::deep_copy(ownedView, particles);
    _ownedParticlesAoS = ownedView;

    _converter.loadDataLayout<Kokkos::View<Particle_T*, HostSpace>, KokkosSoAType<HostSpace>, Particle_T>(particles, ownedSoA, particles.extent(0));
    _converter.storeDataLayout<KokkosSoAType<HostSpace>, Kokkos::View<Particle_T*, HostSpace>, Particle_T>(ownedSoA, particles, ownedView.extent(0));
  }

  virtual void setHaloToTraverse(const Kokkos::View<Particle_T*, HostSpace> &particles) {
    auto haloView = Kokkos::create_mirror_view(DeviceSpace{}, particles);
    Kokkos::deep_copy(haloView, particles);
    _haloParticlesAoS = haloView;
  };

protected:

  Kokkos::View<Particle_T*, DeviceSpace> _ownedParticlesAoS;
  Kokkos::View<Particle_T*, DeviceSpace> _haloParticlesAoS;

  template <typename Space>
  using KokkosSoAType = typename Particle_T::KokkosSoAArraysType<Space>;

  KokkosSoAType<DeviceSpace> _ownedParticleSoA;

private:
  DataLayoutOption _dataLayout;
  utils::KokkosDataLayoutConverter _converter;
};
}