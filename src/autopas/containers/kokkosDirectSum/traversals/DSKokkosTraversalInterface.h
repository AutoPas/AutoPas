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

  void setOwnedToTraverse(utils::KokkosStorage<HostSpace, Particle_T>& input) {

    auto storageLayout = input.getLayout();
    size_t N = input.size();

    if (storageLayout == DataLayoutOption::aos) {
      _ownedParticles.getAoS().resize(N);
      Kokkos::deep_copy(_ownedParticles.getAoS().getView(), input.getAoS().getView());
    }
    else if (storageLayout == DataLayoutOption::soa) {
      _ownedParticles.getSoA().resize(N);

      constexpr size_t tupleSize = input.tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();
      _ownedParticles.getSoA().copyFrom(input.getSoA(), I);
    }
  }

  void setHaloToTraverse(utils::KokkosStorage<HostSpace, Particle_T>& input) {

    auto storageLayout = input.getLayout();
    size_t N = input.size();

    if (storageLayout == DataLayoutOption::aos) {
      _haloParticles.getAoS().resize(N);
      Kokkos::deep_copy(_haloParticles.getAoS().getView(), input.getAoS().getView());
    }
    else if (storageLayout == DataLayoutOption::soa) {
      _haloParticles.getSoA().resize(N);

      constexpr size_t tupleSize = input.tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();
      _haloParticles.getSoA().copyFrom(input.getSoA(), I);
    }
  }

  void retrieveOwned(utils::KokkosStorage<HostSpace, Particle_T>& output) {
    auto storageLayout = _ownedParticles.getLayout();
    size_t N = _ownedParticles.size();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS().resize(N);
      Kokkos::deep_copy(output.getAoS().getView(), _ownedParticles.getAoS().getView());
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA().resize(N);
      constexpr size_t tupleSize = output.tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();
      output.getSoA().copyFrom(_ownedParticles.getSoA(), I);
    }
  }

  void retrieveHalo(utils::KokkosStorage<HostSpace, Particle_T>& output) {
    auto storageLayout = _haloParticles.getLayout();
    size_t N = _haloParticles.size();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS().resize(N);
      Kokkos::deep_copy(output.getAoS().getView(), _haloParticles.getAoS().getView());
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA().resize(N);
      constexpr size_t tupleSize = output.tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();
      output.getSoA().copyFrom(_haloParticles.getSoA(), I);
    }
  }

protected:

  utils::KokkosStorage<DeviceSpace, Particle_T> _ownedParticles;
  utils::KokkosStorage<DeviceSpace, Particle_T> _haloParticles;
};
}