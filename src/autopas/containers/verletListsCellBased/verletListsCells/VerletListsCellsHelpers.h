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
class VerletListsCellsHelpers {
 public:
  /**
   * Cell wise verlet lists: For every cell, a vector of pairs. Each pair maps a particle to a vector of its neighbors.
   */
  using NeighborListsType = std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>;

  /**
   * Pairwise verlet lists: For every cell a vector, for every neighboring cell a vector of particle-neighborlist pairs.
   * Each pair maps a particle to a vector of its neighbor particles.
   * Cells<NeighboringCells<Particle,NeighborParticles>>
   */
  using PairwiseNeighborListsType =
      std::vector<std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>>;

  /**
   * Indicates which type of neighbor list the VLC container is using.
   * To be passed to the VLC traversals in TraversalSelector.
   */
  class VLCTypeOfList {
   public:
    /**
     * Enum value indicating he type of neighbor list.
     */
    enum Value {
      vlc,
      vlp,
    };
  };

};  // class VerletListsCellsHelpers
}  // namespace autopas
