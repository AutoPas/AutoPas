/**
* @file DSKokkosTraversalInterface.h
 * @date 05.11.2025
 * @author Luis Gall
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>
#include "autopas/utilsKokkos/KokkosAoS.h"
#include "autopas/utilsKokkos/KokkosDataLayoutConverter.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {
template <class Particle_T>
class DSKokkosTraversalInterface : public TraversalInterface {

public:

  explicit DSKokkosTraversalInterface(DataLayoutOption dataLayout, bool useNewton3) : TraversalInterface(dataLayout, useNewton3) {};

  virtual ~DSKokkosTraversalInterface() = default;

  void traverseParticles() final {

    performTraversal(_ownedParticles, _ownedParticles);
    performTraversal(_ownedParticles, _haloParticles);
  }

  void setOwnedToTraverse(utilsKokkos::KokkosStorage<Particle_T>& input, DataLayoutOption targetLayout) {

    // No actual copy must be done as we do not want the content to be copied
    _ownedParticles.setIntendedLayout(targetLayout);
    _ownedParticles.setActiveLayout(targetLayout);
    _ownedParticles.overrideSize(input.size());
    _ownedParticles.overrideCapacity(input.getCapacity());

    if (targetLayout == DataLayoutOption::aos) {
      _ownedParticles.getAoS() = input.getAoS();
    }
    else if (targetLayout == DataLayoutOption::soa) {
      _ownedParticles.getSoA() = input.getSoA();
    }
  }

  void setHaloToTraverse(utilsKokkos::KokkosStorage<Particle_T>& input, DataLayoutOption targetLayout) {

    _haloParticles.setIntendedLayout(targetLayout);
    _haloParticles.setActiveLayout(targetLayout);
    _haloParticles.overrideSize(input.size());
    _haloParticles.overrideCapacity(input.getCapacity());

    if (targetLayout == DataLayoutOption::aos) {
      _haloParticles.getAoS() = input.getAoS();
    }
    else if (targetLayout == DataLayoutOption::soa) {
      _haloParticles.getSoA() = input.getSoA();
    }
  }

  void retrieveOwned(utilsKokkos::KokkosStorage<Particle_T>& output) {
    auto storageLayout = _ownedParticles.getActiveLayout();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS() = _ownedParticles.getAoS();
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA() = _ownedParticles.getSoA();
    }
  }

  void retrieveHalo(utilsKokkos::KokkosStorage<Particle_T>& output) {
    auto storageLayout = _haloParticles.getActiveLayout();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS() = _haloParticles.getAoS();
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA() = _haloParticles.getSoA();
    }
  }

#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
  constexpr static bool useHostView = false;
#else
  using DeviceSpace = Kokkos::HostSpace;
  constexpr static bool useHostView = true;
#endif

  using HostSpace = Kokkos::HostSpace;

protected:

  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  virtual void performTraversal(const utilsKokkos::KokkosStorage<Particle_T>& storageA, const utilsKokkos::KokkosStorage<Particle_T>& storageB) = 0;

  utilsKokkos::KokkosStorage<Particle_T> _ownedParticles;
  utilsKokkos::KokkosStorage<Particle_T> _haloParticles;

};
}

#endif