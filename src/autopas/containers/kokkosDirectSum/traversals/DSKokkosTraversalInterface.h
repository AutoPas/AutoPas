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

namespace autopas {
template <class Particle_T>
class DSKokkosTraversalInterface {

public:

  explicit DSKokkosTraversalInterface() {};

  virtual ~DSKokkosTraversalInterface() = default;

  void setOwnedToTraverse(utils::KokkosStorage<Particle_T>& input) {

    auto storageLayout = input.getLayout();
    _ownedParticles.setLayout(storageLayout);
    size_t N = input.size();

    if (storageLayout == DataLayoutOption::aos) {
      _ownedParticles.getAoS() = input.getAoS();
      // Kokkos::deep_copy(_ownedParticles.getAoS().getView(), input.getAoS().getView()); TODO: sync particles
    }
    else if (storageLayout == DataLayoutOption::soa) {
      _ownedParticles.getSoA() = input.getSoA();

      // constexpr size_t tupleSize = input.tupleSize();
      // constexpr auto I = std::make_index_sequence<tupleSize>();
      // _ownedParticles.getSoA().copyFrom(input.getSoA(), I);
    }
  }

  void setHaloToTraverse(utils::KokkosStorage<Particle_T>& input) {

    auto storageLayout = input.getLayout();
    _haloParticles.setLayout(storageLayout);

    if (storageLayout == DataLayoutOption::aos) {
      // _haloParticles.getAoS().resize(N);
      // Kokkos::deep_copy(_haloParticles.getAoS().getView(), input.getAoS().getView()); TODO: sync particles
    }
    else if (storageLayout == DataLayoutOption::soa) {
      // _haloParticles.getSoA().resize(N);
      // constexpr size_t tupleSize = input.tupleSize();
      // constexpr auto I = std::make_index_sequence<tupleSize>();
      // _haloParticles.getSoA().copyFrom(input.getSoA(), I); TODO: sync particles
    }
  }

  void retrieveOwned(utils::KokkosStorage<Particle_T>& output) {
    auto storageLayout = _ownedParticles.getLayout();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS() = _ownedParticles.getAoS();
      // Kokkos::deep_copy(output.getAoS().getView(), _ownedParticles.getAoS().getView()); TODO: sync particles
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA() = _ownedParticles.getSoA();
      // constexpr size_t tupleSize = output.tupleSize();
      // constexpr auto I = std::make_index_sequence<tupleSize>();
      // output.getSoA().copyFrom(_ownedParticles.getSoA(), I); TODO: sync particles
    }
  }

  void retrieveHalo(utils::KokkosStorage<Particle_T>& output) {
    auto storageLayout = _haloParticles.getLayout();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS() = _haloParticles.getAoS();
      // Kokkos::deep_copy(output.getAoS().getView(), _haloParticles.getAoS().getView()); TODO: sync particles
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA() = _haloParticles.getSoA();
      // constexpr size_t tupleSize = output.tupleSize();
      // constexpr auto I = std::make_index_sequence<tupleSize>();
      // output.getSoA().copyFrom(_haloParticles.getSoA(), I); TODO: sync particles
    }
  }

protected:

  utils::KokkosStorage<Particle_T> _ownedParticles;
  utils::KokkosStorage<Particle_T> _haloParticles;
};
}