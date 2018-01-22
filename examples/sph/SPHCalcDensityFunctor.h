//
// Created by seckler on 19.01.18.
//

#ifndef AUTOPAS_SPHCALCDENSITYFUNCTOR_H
#define AUTOPAS_SPHCALCDENSITYFUNCTOR_H

#include "autopas.h"
#include "SPHParticle.h"


class SPHCalcDensityFunctor : public autopas::Functor<SPHParticle> {
public:

    void AoSFunctor(SPHParticle & i, SPHParticle & j) override ;

};


#endif //AUTOPAS_SPHCALCDENSITYFUNCTOR_H
