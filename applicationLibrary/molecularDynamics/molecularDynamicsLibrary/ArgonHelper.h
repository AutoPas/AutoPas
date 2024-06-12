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
 * @tparam T floating point type
 * @param x variable
 * @param order order of the Legendre polynomial
 * @return Lagrange Polynomial of given order
 */
template <class T>
[[nodiscard]] constexpr T Legendre(const T &x, const size_t order) {
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
 * @tparam T floating point type
 * @param x variable
 * @param order order of the Legendre polynomial
 * @return derivative w.r.t. x of the Legendre polynomial
 */
template <class T>
[[nodiscard]] constexpr T derivativeLegendre(const T& x, size_t order) {
  auto Legendre_order_minus_1 = (order == 0) ? 0 : Legendre(x, order - 1);
  auto num{Legendre_order_minus_1 - x * Legendre(x, order)};
  return order * num / (1 - Math::pow<2>(x));
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

template <size_t wrt>
[[nodiscard]] std::array<double, 3> P(const size_t A, const size_t B, const size_t C, const Cosine& cosI, const Cosine& cosJ, const Cosine& cosK) {
  auto firstTerm = derivativeLegendre(cosI.getCos(), A) * cosI.derive_wrt<wrt>() * Legendre(cosJ.getCos(), B) * Legendre(cosK.getCos(), C);
  auto secondTerm = Legendre(cosI.getCos(), A) * derivativeLegendre(cosJ.getCos(), B) * cosJ.derive_wrt<wrt>() * Legendre(cosK.getCos(), C);
  auto thirdTerm = Legendre(cosI.getCos(), A) * Legendre(cosJ.getCos(), B) * derivativeLegendre(cosK.getCos(), C) * cosK.derive_wrt<wrt>();
  return firstTerm + secondTerm + thirdTerm;
}

template <size_t wrt>
[[nodiscard]] std::array<double, 3> Permutation_deriv_wrt(const size_t A, const size_t B, const size_t C, const Cosine& cosI, const Cosine& cosJ, const Cosine& cosK) {
  return P<wrt>(A, B, C, cosI, cosJ, cosK) + P<wrt>(A, C, B, cosI, cosJ, cosK)  + P<wrt>(B, A, C, cosI, cosJ, cosK)
      + P<wrt>(B, C, A, cosI, cosJ, cosK) + P<wrt>(C, A, B, cosI, cosJ, cosK) + P<wrt>(C, B, A, cosI, cosJ, cosK);
}

} // namespace autopas::utils::ArgonMath
