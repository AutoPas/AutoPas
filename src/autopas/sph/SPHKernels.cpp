/**
 * @file SPHKernels.cpp
 * @author seckler
 * @date 22.01.18
 */
#include "autopas/sph/SPHKernels.h"

#include <algorithm>
#include <cmath>

uint64_t autopas::sph::SPHKernels::getFlopsW() {
  uint64_t flops = 0;
  flops += 1;      // calculating H
  flops += 5;      // dot product for s
  flops += 1 + 1;  // s (except dot product)
  flops += 1;      // calculating s1 and s2 (either for s1 or s2 one flop will be
                   // necessary
  flops += 6 + 4;  // calculating r_value
  return flops;
}
