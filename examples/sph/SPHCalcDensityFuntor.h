//
// Created by seckler on 19.01.18.
//

#ifndef AUTOPAS_SPHCALCDENSITYFUNTOR_H
#define AUTOPAS_SPHCALCDENSITYFUNTOR_H

#include "autopas.h"
#include "SPHParticle.h"
//#include "particles/MoleculeLJ.h"
#include <array>


class SPHCalcDensityFunctor : public autopas::Functor<SPHParticle> {
public:


    void AoSFunctor(SPHParticle & i, SPHParticle & j) override ;
//    void operator()(const EP *const ep_i, const PS::S32 Nip,
//                    const EP *const ep_j, const PS::S32 Njp,
//                    Dens *const dens) {
//        for (PS::S32 i = 0; i < Nip; ++i) {
//            dens[i].clear();
//            for (PS::S32 j = 0; j < Njp; ++j) {
//                const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
//                dens[i].dens += ep_j[j].mass * W(dr, ep_i[i].smth);
//            }
//        }
//    }

};


#endif //AUTOPAS_SPHCALCDENSITYFUNTOR_H
