/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 11/06/24
 */

#pragma once

#include "CosineHandle.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::utils::ArrayMath::Argon {

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

}   // namespace autopas::utils::ArrayMath::Argon
