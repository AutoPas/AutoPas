/**
 * @file KokkosParticleCell.h
 * @date 20.10.2021
 * @author lgaertner
 */

#pragma once

#include <Kokkos_Core.hpp>
#include <array>

#include "autopas/cells/ParticleCell.h"

namespace autopas {

/**
 * Dummy particle cell to use CellBlock3D with Kokkos Containers
 * @tparam Particle
 */
template <class Particle>
class KokkosParticleCell {
 public:
  /**
   * The particle type for this cell.
   */
  using ParticleType = Particle;

  KokkosParticleCell() : particlesPtr(nullptr), begin(0ul), cellSize(0ul){};

  std::array<size_t, 2> getRange() { return {begin, begin + cellSize}; }

  size_t begin;
  size_t cellSize;
  Kokkos::View<Particle *> *particlesPtr;
};
}  // namespace autopas