//
// Created by seckler on 19.01.18.
//

#include "SPHCalcDensityFunctor.h"
#include "SPHKernels.h"
#include "utils/arrayMath.h"

using namespace autopas::sph;

unsigned long SPHCalcDensityFunctor::getNumFlopsPerKernelCall() {
  unsigned long flops = 0;
  flops += 3;                            // calculating dr
  flops += 2 * SPHKernels::getFlopsW();  // flops for calling W
  flops += 2 * 1;                        // calculating density
  flops += 2 * 1;                        // adding density
  return flops;
}