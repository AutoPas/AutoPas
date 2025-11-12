/**
* @file DSKokkosTraversalInterface.h
 * @date 05.11.2025
 * @author Luis Gall
 */

#pragma once

#include <Kokkos_Core.hpp>

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

  explicit DSKokkosTraversalInterface() {};

  virtual ~DSKokkosTraversalInterface() = default;

  // TODO: data Copy
  virtual void setOwnedToTraverse(const Kokkos::View<Particle_T*, HostSpace> &particles) {
    auto ownedView = Kokkos::create_mirror_view(DeviceSpace{}, particles);
    Kokkos::deep_copy(ownedView, particles);
    _ownedParticles = ownedView;
  };
  virtual void setHaloToTraverse(const Kokkos::View<Particle_T*, HostSpace> &particles) {
    auto haloView = Kokkos::create_mirror_view(DeviceSpace{}, particles);
    Kokkos::deep_copy(haloView, particles);
    _haloParticles = haloView;
  };

protected:

  Kokkos::View<Particle_T*, DeviceSpace> _ownedParticles;
  Kokkos::View<Particle_T*, DeviceSpace> _haloParticles;
};
}