//
// Created by seckler on 19.01.18.
//

#include "SPHCalcDensityFunctor.h"
#include "utils/arrayMath.h"
#include "SPHKernels.h"

using namespace autopas::sph;

void SPHCalcDensityFunctor::AoSFunctor(SPHParticle &i, SPHParticle &j) {
    const std::array<double, 3> dr = arrayMath::sub(j.getR(), i.getR());  // ep_j[j].pos - ep_i[i].pos;
    const double density = j.getMass() * W(dr, i.getSmth());  // ep_j[j].mass * W(dr, ep_i[i].smth)
    i.addDensity(density);
}