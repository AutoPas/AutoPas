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
  explicit CosineHandle(const DisplacementHandle &displacementAB, const DisplacementHandle &displacementAC);

  [[nodiscard]] double getCos() const { return cos_; }

  template <size_t wrt>
  [[nodiscard]] nabla derive_wrt();

 private:
  DisplacementHandle displacementAB_;
  DisplacementHandle displacementAC_;
  std::array<double, 3> AB_;
  std::array<double, 3> AC_;
  double cos_;
  size_t id_;
};

}  // namespace autopas::utils::ArrayMath::Argon