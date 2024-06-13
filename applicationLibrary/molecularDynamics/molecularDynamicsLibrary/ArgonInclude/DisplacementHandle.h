/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 13/06/24
 */

#pragma once

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::utils::ArrayMath::Argon {
using nabla = std::array<double, 3>;

class DisplacementHandle {
 public:
  DisplacementHandle() = default;
  explicit DisplacementHandle(const std::array<double, 3>& positionStartVertex, const std::array<double, 3>& positionEndVertex,
                        const size_t& idStartVertex, const size_t& idEndVertex);

  [[nodiscard]] size_t getIdStartVertex() const { return idStartVertex_; }
  [[nodiscard]] size_t getIdEndVertex() const { return idEndVertex_; }
  [[nodiscard]] std::array<double, 3> getDisplacement() const { return  displacement_; }
  [[nodiscard]] DisplacementHandle getInv() const;

  template<size_t wrt>
  [[nodiscard]] nabla derive_wrt();

 private:
  std::array<double, 3> positionStartVertex_;
  std::array<double, 3> positionEndVertex_;
  std::array<double, 3> displacement_;
  size_t idStartVertex_;
  size_t idEndVertex_;
};

}   // namespace autopas::utils::ArrayMath::Argon