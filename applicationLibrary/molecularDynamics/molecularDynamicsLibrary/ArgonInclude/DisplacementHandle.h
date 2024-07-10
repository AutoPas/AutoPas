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

  /**
   * Constructor for DisplacementHandle. It stores the ids and the positions of the start and end vertices, as well as
   * the displacement vector between the two (i.e. positionEndVertex - positionStartVertex)
   * @param positionStartVertex position of the start vertex
   * @param positionEndVertex position of the end vertex
   * @param idStartVertex id of the start vertex
   * @param idEndVertex id of the end vertex
   */
  explicit DisplacementHandle(const std::array<double, 3> &positionStartVertex,
                              const std::array<double, 3> &positionEndVertex, const size_t &idStartVertex,
                              const size_t &idEndVertex);

  /**
   *
   * @return id of the start vertex
   */
  [[nodiscard]] size_t getIdStartVertex() const { return idStartVertex_; }

  /**
   *
   * @return id of the end vertex
   */
  [[nodiscard]] size_t getIdEndVertex() const { return idEndVertex_; }

  /**
   *
   * @return array obtained from positionEndVertex - positionStartVertex
   */
  [[nodiscard]] std::array<double, 3> getDisplacement() const { return displacement_; }

  /**
   *
   * @return DisplacementHandle with inverted start and end vertices
   */
  [[nodiscard]] DisplacementHandle getInv() const;

  /**
   *
   * @tparam ID id of the particle with respect to which we are computing the derivative
   * @return derivative of the cosine displacement_ w.r.t. ID
   */
  template <size_t ID>
  [[nodiscard]] nabla derive_wrt() const;

 private:
  std::array<double, 3> positionStartVertex_;
  std::array<double, 3> positionEndVertex_;
  std::array<double, 3> displacement_;
  size_t idStartVertex_;
  size_t idEndVertex_;
};

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

template <size_t ID>
[[nodiscard]] nabla DisplacementHandle::derive_wrt() const {
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