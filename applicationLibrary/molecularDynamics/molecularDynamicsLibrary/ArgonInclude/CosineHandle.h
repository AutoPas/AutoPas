/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 13/06/24
 */

#pragma once

#include "DisplacementHandle.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::utils::ArrayMath::Argon {

class CosineHandle {
 public:
  /**
   * constructor for the CosineHandle. Constructs the CosineHandle if displacementAB.idStartVertex ==
   * displacementAC.idStartVertex
   * @param displacementAB
   * @param displacementAC
   */
  explicit CosineHandle(const DisplacementHandle &displacementAB, const DisplacementHandle &displacementAC);

  /**
   *
   * @return cosine of the angle between displacementAB.displacement_ and displacementAC.displacement_
   */
  [[nodiscard]] double getCos() const { return cos_; }

  /**
   *
   * @param ID id of the particle with respect to which we are computing the derivative
   * @return derivative of the cosine cos_ w.r.t. ID
   */
  [[nodiscard]] nabla derive_wrt(const size_t ID) const;

 private:
  DisplacementHandle displacementAB_;
  DisplacementHandle displacementAC_;
  std::array<double, 3> AB_;
  std::array<double, 3> AC_;
  double cos_;
  size_t id_;
};

}  // namespace autopas::utils::ArrayMath::Argon