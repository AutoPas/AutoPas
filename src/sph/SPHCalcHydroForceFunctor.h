//
// Created by seckler on 22.01.18.
//

#ifndef AUTOPAS_SPHCALCHYDROFORCEFUNCTOR_H
#define AUTOPAS_SPHCALCHYDROFORCEFUNCTOR_H

#include "SPHParticle.h"
#include "autopas.h"

namespace autopas {
namespace sph {
class SPHCalcHydroForceFunctor : public autopas::Functor<SPHParticle> {
 public:
  void AoSFunctor(SPHParticle &i, SPHParticle &j) override;

  static unsigned long getNumFlopsPerKernelCall();
};
}  // namespace sph
}  // namespace autopas
#endif  // AUTOPAS_SPHCALCHYDROFORCEFUNCTOR_H
