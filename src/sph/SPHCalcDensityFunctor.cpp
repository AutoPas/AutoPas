//
// Created by seckler on 19.01.18.
//

#include "SPHCalcDensityFunctor.h"
#include "utils/arrayMath.h"
#include "SPHKernels.h"

using namespace autopas::sph;

void SPHCalcDensityFunctor::AoSFunctor(SPHParticle &i, SPHParticle &j) {
    const std::array<double, 3> dr = arrayMath::sub(j.getR(), i.getR());  // ep_j[j].pos - ep_i[i].pos;
    const double density = j.getMass() * SPHKernels::W(dr, i.getSmth());  // ep_j[j].mass * W(dr, ep_i[i].smth)
    i.addDensity(density);
}

unsigned long SPHCalcDensityFunctor::getNumFlopsPerKernelCall(){
    unsigned long flops = 0;
    flops += 3;  // calculating dr
    flops += SPHKernels::getFlopsW();  // flops for calling W
    flops += 1;  // calculating density
    flops += 1;  // adding density
    return flops;
}