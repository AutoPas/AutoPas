/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 13/06/24
 */

#pragma once

#include "DisplacementHandle.h"

namespace autopas::utils::ArrayMath::Argon {

DisplacementHandle::DisplacementHandle(const std::array<double, 3> &positionStartVertex,
                                       const std::array<double, 3> &positionEndVertex, const size_t &idStartVertex,
                                       const size_t &idEndVertex)
    : positionStartVertex_(positionStartVertex),
      positionEndVertex_(positionEndVertex),
      displacement_(positionEndVertex - positionStartVertex),
      idStartVertex_(idStartVertex),
      idEndVertex_(idEndVertex) {}

[[nodiscard]] DisplacementHandle DisplacementHandle::getInv() const {
  return DisplacementHandle(positionEndVertex_, positionStartVertex_, idEndVertex_, idStartVertex_);
}

//template <size_t wrt>
[[nodiscard]] nabla DisplacementHandle::derive_wrt(size_t ID) const {
  auto moduloDisplacement{L2Norm(displacement_)};
  if (ID != idStartVertex_ && ID != idEndVertex_) {
    return std::array<double, 3>{{0, 0, 0}};
  } else if (ID == idStartVertex_) {
    return std::array<double, 3>{0, 0, 0} - displacement_ / moduloDisplacement;
  } else if (ID == idEndVertex_) {
    return displacement_ / moduloDisplacement;
  }
}
}  // namespace autopas::utils::ArrayMath::Argon
