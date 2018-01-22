//
// Created by seckler on 22.01.18.
//

#ifndef AUTOPAS_SPHCALCHYDROFORCEFUNCTOR_H
#define AUTOPAS_SPHCALCHYDROFORCEFUNCTOR_H

#include "autopas.h"
#include "SPHParticle.h"

namespace autopas {
    namespace sph {
        class SPHCalcHydroForceFunctor : public autopas::Functor<SPHParticle> {

            void AoSFunctor(SPHParticle &i, SPHParticle &j) override;

        };
    }  // namespace autopas
}  // namespace autopas
#endif //AUTOPAS_SPHCALCHYDROFORCEFUNCTOR_H
