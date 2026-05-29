/**
 * @file VerletListsKokkosTraversalInterface.h
 * @date 28.05.2026
 * @author Franziska Duhr
 * @note This file was created from DSKokkosTraversalInterface.h. May want to refactor to avoid code duplication in the future.
 */

#pragma once

#include <Kokkos_Core.hpp>
#include "autopas/utils/KokkosAoS.h"
#include "autopas/utils/KokkosDataLayoutConverter.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {
template <class Particle_T>
class VerletListsKokkosTraversalInterface {

public:

  explicit VerletListsKokkosTraversalInterface() {};

  virtual ~VerletListsKokkosTraversalInterface() = default;

  void setOwnedToTraverse(utils::KokkosStorage<Particle_T>& input) {

    auto storageLayout = input.getLayout();
    _ownedParticles.setLayout(storageLayout);

    if (storageLayout == DataLayoutOption::aos) {
      _ownedParticles.getAoS() = input.getAoS();
    }
    else if (storageLayout == DataLayoutOption::soa) {
      _ownedParticles.getSoA() = input.getSoA();
    }
  }

  void setHaloToTraverse(utils::KokkosStorage<Particle_T>& input) {

    auto storageLayout = input.getLayout();
    _haloParticles.setLayout(storageLayout);

    if (storageLayout == DataLayoutOption::aos) {
      _haloParticles.getAoS() = input.getAoS();
    }
    else if (storageLayout == DataLayoutOption::soa) {
      _haloParticles.getSoA() = input.getSoA();
    }
  }

  void retrieveOwned(utils::KokkosStorage<Particle_T>& output) {
    auto storageLayout = _ownedParticles.getLayout();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS() = _ownedParticles.getAoS();
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA() = _ownedParticles.getSoA();
    }
  }

  void retrieveHalo(utils::KokkosStorage<Particle_T>& output) {
    auto storageLayout = _haloParticles.getLayout();

    if (storageLayout == DataLayoutOption::aos) {
      output.getAoS() = _haloParticles.getAoS();
    }
    else if (storageLayout == DataLayoutOption::soa) {
      output.getSoA() = _haloParticles.getSoA();
    }
  }

  void setNeighborList(Kokkos::View<size_t*> offsets, Kokkos::View<size_t*> entries) {
    _neighborListOffsets = offsets;
    _neighborListEntries = entries;
  }

protected:

  utils::KokkosStorage<Particle_T> _ownedParticles;
  utils::KokkosStorage<Particle_T> _haloParticles;

  Kokkos::View<size_t*> _neighborListOffsets;
  Kokkos::View<size_t*> _neighborListEntries;
};
}
