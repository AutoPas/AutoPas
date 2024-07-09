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

[[nodiscard]] nabla CosineHandle::derive_wrt(const size_t ID) const {
  if (id_ != ID && displacementAB_.getIdEndVertex() != ID && displacementAC_.getIdEndVertex() != ID) {
    return std::array<double, 3>{{0, 0, 0}};
  } else if (ID == id_) {
    auto firstTerm{cos_ / ArrayMath::dot(AB_, AB_) - 1. / ArrayMath::dot(AB_, AC_)};
    auto secondTerm{cos_ / ArrayMath::dot(AC_, AC_) - 1. / ArrayMath::dot(AB_, AC_)};
    return AB_ * firstTerm + AC_ * secondTerm;
  } else if (ID == displacementAB_.getIdEndVertex()) {
    auto firstTerm{-cos_ / ArrayMath::dot(AB_, AB_)};
    auto secondTerm{1. / ArrayMath::dot(AB_, AC_)};
    return AB_ * firstTerm + AC_ * secondTerm;
  } else if (ID == displacementAC_.getIdEndVertex()) {
    auto firstTerm{-cos_ / ArrayMath::dot(AC_, AC_)};
    auto secondTerm{1. / ArrayMath::dot(AB_, AC_)};
    return AC_ * firstTerm + AB_ * secondTerm;
  }
}

}  // namespace autopas::utils::ArrayMath::Argon