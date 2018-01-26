//
// Created by seckler on 19.01.18.
//

#include "SPHCalcDensityFunctor.h"
#include "utils/arrayMath.h"
#include "SPHKernels.h"

using namespace autopas::sph;


unsigned long SPHCalcDensityFunctor::getNumFlopsPerKernelCall(){
    unsigned long flops = 0;
    flops += 3;  // calculating dr
    flops += SPHKernels::getFlopsW();  // flops for calling W
    flops += 1;  // calculating density
    flops += 1;  // adding density
    return flops;
}