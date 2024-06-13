/**
* @file ArgonFunctor.h
* @author I. Angelucci
* @date 13/06/24
*/

#pragma once

#include "CosineHandle.h"

namespace autopas::utils::ArrayMath::Argon {

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

  /**
   *
   * @tparam wrt particle with respect to which we are computing the derivative
   * @return Derivative of cos_ with respect to particle wrt
   */
  template<size_t wrt>
  [[nodiscard]] nabla CosineHandle::derive_wrt() {
    if(id_ != wrt && displacementAB_.getIdEndVertex() != wrt && displacementAC_.getIdEndVertex() != wrt) {
      return std::array<double, 3>{{0, 0, 0}};
    }
    else if (wrt == id_) {
      auto firstTerm{cos_ / ArrayMath::dot(AB_, AB_) - 1. / ArrayMath::dot(AB_, AC_)};
      auto secondTerm{cos_ / ArrayMath::dot(AC_, AC_) - 1. / ArrayMath::dot(AB_, AC_)};
      return AB_ * firstTerm + AC_ * secondTerm;
    }
    else if (wrt == displacementAB_.getIdEndVertex()) {
      auto firstTerm{-cos_ / ArrayMath::dot(AB_, AB_)};
      auto secondTerm{1./ ArrayMath::dot(AB_, AC_)};
      return AB_ * firstTerm + AC_ * secondTerm;
    }
    else if (wrt == displacementAC_.getIdEndVertex()) {
      auto firstTerm{-cos_ / ArrayMath::dot(AC_, AC_)};
      auto secondTerm{1./ ArrayMath::dot(AB_, AC_)};
      return AC_ * firstTerm + AB_ * secondTerm;
    }
  }

}   // namespace autopas::utils::ArrayMath::Argon