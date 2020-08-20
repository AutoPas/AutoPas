/**
 * @file VerletListTypeDefinitions.h
 * @author seckler
 * @date 13.07.20
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/utils/SoAType.h"

namespace autopas {
/**
 * Class of helpers for verlet lists. This is a class and not a namespace to allow templating.
 * @tparam Particle
 */
template <class Particle>
class VerletListTypeDefinitions {
 public:
  /// using declaration for soa's of verlet list's linked cells (only id and position needs to be stored)
  using PositionSoAArraysType = typename utils::SoAType<Particle *, double, double, double>::Type;

  /// using declaration for verlet-list particle cell type
  using VerletListParticleCellType = FullParticleCell<Particle, PositionSoAArraysType>;
};
}  // namespace autopas
