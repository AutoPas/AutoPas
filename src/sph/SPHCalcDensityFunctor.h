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
class SPHCalcDensityFunctor : public Functor<SPHParticle> {
 public:
  inline void AoSFunctor(SPHParticle &i, SPHParticle &j) override {
    const std::array<double, 3> dr =
        arrayMath::sub(j.getR(), i.getR());  // ep_j[j].pos - ep_i[i].pos;
    const double density =
        j.getMass() *
        SPHKernels::W(dr, i.getSmth());  // ep_j[j].mass * W(dr, ep_i[i].smth)
    i.addDensity(density);
  }
  static unsigned long getNumFlopsPerKernelCall();
};
}  // namespace sph
}  // namespace autopas

#endif  // AUTOPAS_SPHCALCDENSITYFUNCTOR_H
