/**
 * @file NewVerletListHelpers.h
 * @author Luis Gall
 * @date 04.05.2023
 *
 * oriented on
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

template <class Particle>
class NewVerletListHelpers {

 public:

  /**'
   * Static Neighbor Lists
   */
  using StaticNeighborListsType = std::unordered_map<Particle *, std::vector<Particle *>>;

  /**
   * Dynamic Neighbor Lists
   */
  using DynamicNeighborListsType = std::unordered_map<Particle *, std::pair<std::vector<Particle *>, std::array<double, 3>>>;

  class VLBuildType {
   public:
    enum Value {
      aosBuild,
      soaBuild
    };
  };
};

}