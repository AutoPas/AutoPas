//
// Created by seckler on 22.01.18.
//
#pragma once

#include "SPHParticle.h"
#include "autopasIncludes.h"

namespace autopas {
namespace sph {
/**
 * Class that defines the hydrodynamic force functor.
 * It is used to calculate the force based on the given SPH kernels.
 */
class SPHCalcHydroForceFunctor
    : public autopas::Functor<
          SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> {
 public:
  /**
   * Calculates the contribution of the interaction of particle i and j to the
   * hydrodynamic force.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  void AoSFunctor(SPHParticle &i, SPHParticle &j, bool newton3 = true) override;

  /**
   * Get the number of floating point operations used in one full kernel call
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall();
};
}  // namespace sph
}  // namespace autopas

