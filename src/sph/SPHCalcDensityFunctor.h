//
// Created by seckler on 19.01.18.
//

#ifndef AUTOPAS_SPHCALCDENSITYFUNCTOR_H
#define AUTOPAS_SPHCALCDENSITYFUNCTOR_H

#include "SPHKernels.h"
#include "SPHParticle.h"
#include "autopas.h"

namespace autopas {
namespace sph {
/**
 * Class that defines the density functor.
 * It is used to calculate the density based on the given SPH kernel.
 */
class SPHCalcDensityFunctor : public Functor<SPHParticle> {
 public:
  inline void AoSFunctor(SPHParticle &i, SPHParticle &j) override {
    const std::array<double, 3> dr =
        arrayMath::sub(j.getR(), i.getR());  // ep_j[j].pos - ep_i[i].pos;
    const double density =
        j.getMass() *
            SPHKernels::W(dr, i.getSmoothingLength());  // ep_j[j].mass * W(dr, ep_i[i].smth)
    i.addDensity(density);
  }
  /**
   * Get the number of floating point operations used in one full kernel call
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall();
};
}  // namespace sph
}  // namespace autopas

#endif  // AUTOPAS_SPHCALCDENSITYFUNCTOR_H
