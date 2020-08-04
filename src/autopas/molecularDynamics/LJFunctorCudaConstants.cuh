/**
 * @file LJFunctorCudaConstants.cuh
 * @author seckler
 * @date 19.06.20
 */

#pragma once

#include "autopas/pairwiseFunctors/FunctorCuda.cuh"

namespace autopas {

/**
 * Stores all constants needed for the calculation
 * @tparam floatType of constants
 */
template <typename floatType>
class LJFunctorConstants : public FunctorCudaConstants<floatType> {
 public:
  LJFunctorConstants() {}
  LJFunctorConstants(floatType csq, floatType ep24, floatType sqs, floatType sh6)
      : cutoffsquare(csq), epsilon24(ep24), sigmasquare(sqs), shift6(sh6) {}
  floatType cutoffsquare;
  floatType epsilon24;
  floatType sigmasquare;
  floatType shift6;
};

}  // namespace autopas