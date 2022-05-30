/**
 * @file VerletListsCellsHelpers.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"

namespace autopas {

/**
 * class of helpers for verlet lists
 * @tparam Particle
 */
template <class Particle>
class PseudoVerletListHelpers {
 public:
  /**
   * Cell wise pseudo verlet lists: For every cell, a vector of pairs. Each pair maps a particle to its furthest potential neighbor.
   */
  using PairwiseNeighborListsType = std::vector<std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>>;

  /**
   * Indicates which build functor should be used for the generation of the neighbor list.
   * To be passed to the generator functors in the neighbor lists.
   */
  class PVLBuildType {
   public:
    /**
     * Enum value indicating whether the AoS functor or the SoA functors will be used to generate
     * the neighbor list.
     */
    enum Value {
      aosBuild,
      soaBuild,
    };
  };

};  // class PseudoVerletListHelpers
}  // namespace autopas
