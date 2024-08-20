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
   * @tparam ID id of the particle with respect to which we are computing the derivative
   * @return derivative of the cosine cos_ w.r.t. ID
   */
  template <size_t ID>
  [[nodiscard]] nabla derive_wrt() const;

 private:
  DisplacementHandle displacementAB_;
  DisplacementHandle displacementAC_;
  std::array<double, 3> AB_;
  std::array<double, 3> AC_;
  double cos_;
  size_t id_;
};

CosineHandle::CosineHandle(const DisplacementHandle &displacementAB, const DisplacementHandle &displacementAC) {
  if (displacementAB.getIdStartVertex() != displacementAC.getIdStartVertex()) {
    throw autopas::utils::ExceptionHandler::AutoPasException("cannot build cosine");
  }
  displacementAB_ = displacementAB;
  displacementAC_ = displacementAC;
  AB_ = displacementAB_.getDisplacement();
  AC_ = displacementAC_.getDisplacement();
  cos_ = ArrayMath::dot(AB_, AC_) / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AC_));
  id_ = displacementAB.getIdStartVertex();
}

template <size_t ID>
[[nodiscard]] nabla CosineHandle::derive_wrt() const {
  if (ID == id_) {
    auto firstTerm{cos_ / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AB_)) -
                   1. / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AC_))};
    auto secondTerm{cos_ / (ArrayMath::L2Norm(AC_) * ArrayMath::L2Norm(AC_)) -
                    1. / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AC_))};
    return AB_ * firstTerm + AC_ * secondTerm;
  } else if (ID == displacementAB_.getIdEndVertex()) {
    auto firstTerm{-cos_ / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AB_))};
    auto secondTerm{1. / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AC_))};
    return AB_ * firstTerm + AC_ * secondTerm;
  } else if (ID == displacementAC_.getIdEndVertex()) {
    auto firstTerm{-cos_ / (ArrayMath::L2Norm(AC_) * ArrayMath::L2Norm(AC_))};
    auto secondTerm{1. / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AC_))};
    return AC_ * firstTerm + AB_ * secondTerm;
  }
  return std::array<double, 3>{{0, 0, 0}};
}

}  // namespace autopas::utils::ArrayMath::Argon