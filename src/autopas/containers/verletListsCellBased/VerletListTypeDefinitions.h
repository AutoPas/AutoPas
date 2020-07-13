/**
 * @file VerletListTypeDefinitions.h
 * @author seckler
 * @date 13.07.20
 */

#pragma once

#include "autopas/utils/SoAType.h"
#include "autopas/cells/FullParticleCell.h"

namespace autopas {
/**
 * class of helpers for verlet lists
 * @tparam Particle
 */
template <class Particle>
class VerletListTypeDefinitions {
 public:
  /// using declaration for soa's of verlet list's linked cells (only id and position needs to be stored)
  using SoAArraysType = typename utils::SoAType<Particle *, double, double, double>::Type;

  /// using declaration for verlet-list particle cell type
  using VerletListParticleCellType = FullParticleCell<Particle, SoAArraysType>;
};
}  // namespace autopas