
#pragma once
#include "autopas/containers/CellBasedParticleContainer.h"

namespace autopas::kokkos {
template <class Particle_T>
class DirectSum : public CellBasedParticleContainer<> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle_T>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle_T;
};
}  // namespace autopas::kokkos