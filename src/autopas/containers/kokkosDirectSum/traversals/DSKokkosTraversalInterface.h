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

protected:

  utils::KokkosStorage<Particle_T> _ownedParticles;
  utils::KokkosStorage<Particle_T> _haloParticles;
};
}