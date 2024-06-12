/**
* @file ArgonFunctor.h
* @author I. Angelucci
* @date 11/06/24
 */

#pragma once

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::utils::ArrayMath::ArgonMath {

/**
 *
 * @tparam order order of the Legendre polynomial
 * @tparam T floating point type
 * @param x variable
 * @return Lagrange Polynomial of given order
 */
template <size_t order, class T>
[[nodiscard]] constexpr T Legendre(const T &x) {
  if (order > 6) {
    throw ExceptionHandler::AutoPasException(
        "Argon Simulation: asked to compute Legendre polynomial of order higher than 6.");
  }
  switch (order) {
    case 0: return 1;
    case 1: return x;
    case 2: return 0.5 * (3 * Math::pow<2>(x) - 1);
    case 3: return 0.5 * (5 * Math::pow<3>(x) - 3 * x);
    case 4: return 0.125 * (35 * Math::pow<4>(x) - 30 * Math::pow<2>(x) + 3);
    case 5: return 0.125 * (63 * Math::pow<5>(x) - 70 * Math::pow<3>(x) + 15 * x);
    case 6: return 0.0625 * (231 * Math::pow<6>(x) - 315 * Math::pow<4>(x) + 105 * Math::pow<2>(x) - 5);
  }
}

/**
 *
 * @tparam order order of the Legendre polynomial
 * @tparam T floating point type
 * @param x variable
 * @return derivative w.r.t. x of the Legendre polynomial
 */
template <size_t order, class T>
[[nodiscard]] constexpr T derivativeLegendre(const T &x) {
  auto Legendre_order_minus_1 = (order == 0) ? 0 : Legendre<order - 1>(x);
  auto num{Legendre_order_minus_1 - Math::safeMul(x, Legendre<order>(x))};
  return Math::safeMul(order, num) / (1 - Math::pow<2>(x));
}

class Displacement {
 public:
  Displacement() = default;
  explicit Displacement(const std::array<double, 3>& positionStartVertex, const std::array<double, 3>& positionEndVertex,
               const size_t& idStartVertex, const size_t& idEndVertex) : positionStartVertex_(positionStartVertex), positionEndVertex_(positionEndVertex),
                                                                       idStartVertex_(idStartVertex), idEndVertex_(idEndVertex) {};

  [[nodiscard]] size_t getIdStartVertex() const { return idStartVertex_; }
  [[nodiscard]] size_t getIdEndVertex() const { return idEndVertex_; }
  [[nodiscard]] std::array<double, 3> getR() const { return  positionEndVertex_ - positionStartVertex_; }
  [[nodiscard]] Displacement getInv() const { return Displacement(positionEndVertex_, positionStartVertex_, idEndVertex_, idStartVertex_); }

 private:
  std::array<double, 3> positionStartVertex_;
  std::array<double, 3> positionEndVertex_;
  size_t idStartVertex_;
  size_t idEndVertex_;
};

class Cosine {
 public:
  explicit Cosine(const Displacement& displacementAB, const Displacement& displacementAC) {
    if (displacementAB.getIdStartVertex() != displacementAC.getIdStartVertex()) {
      throw autopas::utils::ExceptionHandler::AutoPasException("cannot build cosine");
    }
    displacementAB_ = displacementAB;
    displacementAC_ = displacementAC;
    AB_ = displacementAB_.getR();
    AC_ = displacementAC_.getR();
    cos_ = ArrayMath::dot(AB_, AC_) / (ArrayMath::L2Norm(AB_) * ArrayMath::L2Norm(AC_));
    id_ = displacementAB.getIdStartVertex();
  };

  [[nodiscard]] double getCos() const { return cos_; }

  template<size_t wrt>
  [[nodiscard]] std::array<double, 3> derive_wrt() {
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

 private:
  Displacement displacementAB_;
  Displacement displacementAC_;
  std::array<double, 3> AB_;
  std::array<double, 3> AC_;
  double cos_;
  size_t id_;
};

} // namespace autopas::utils::ArgonMath
