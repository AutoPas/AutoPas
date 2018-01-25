//
// Created by seckler on 19.01.18.
//

#ifndef AUTOPAS_SPHCALCDENSITYFUNCTOR_H
#define AUTOPAS_SPHCALCDENSITYFUNCTOR_H

#include "autopas.h"
#include "SPHParticle.h"

namespace autopas {
    namespace sph {
        class SPHCalcDensityFunctor : public Functor<SPHParticle> {
        public:

            void AoSFunctor(SPHParticle &i, SPHParticle &j) override;
            static unsigned long getNumFlopsPerKernelCall();
        };
    }  // namespace autopas
}  // namespace autopas

#endif //AUTOPAS_SPHCALCDENSITYFUNCTOR_H
