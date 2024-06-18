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
    case 0:
      return 1;
    case 1:
      return x;
    case 2:
      return 0.5 * (3 * Math::pow<2>(x) - 1);
    case 3:
      return 0.5 * (5 * Math::pow<3>(x) - 3 * x);
    case 4:
      return 0.125 * (35 * Math::pow<4>(x) - 30 * Math::pow<2>(x) + 3);
    case 5:
      return 0.125 * (63 * Math::pow<5>(x) - 70 * Math::pow<3>(x) + 15 * x);
    case 6:
      return 0.0625 * (231 * Math::pow<6>(x) - 315 * Math::pow<4>(x) + 105 * Math::pow<2>(x) - 5);
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
[[nodiscard]] constexpr T derivativeLegendre(const T &x, size_t order) {
  auto Legendre_order_minus_1 = (order == 0) ? 0 : Legendre(x, order - 1);
  auto num{Legendre_order_minus_1 - x * Legendre(x, order)};
  return order * num / (1 - Math::pow<2>(x));
}

/**
 *
 * @param A order of the Legendre Polynomial for cosI
 * @param B order of the Legendre Polynomial for cosJ
 * @param C order of the Legendre Polynomial for cosK
 * @param cosI
 * @param cosJ
 * @param cosK
 * @return LegendrePolynomial_A(cosI) * LegendrePolynomial_B(cosJ) * LegendrePolynomial_C(cosK)
 */
[[nodiscard]] double P(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                       const CosineHandle &cosJ, const CosineHandle &cosK) {
  return Legendre(cosI.getCos(), A) * Legendre(cosJ.getCos(), B) * Legendre(cosK.getCos(), C);
}

/**
 *
 * @tparam wrt particle with respect to which we are computing the derivative
 * @param A
 * @param B
 * @param C
 * @param cosI
 * @param cosJ
 * @param cosK
 * @return derivative of P (defined above)
 */
template <size_t wrt>
[[nodiscard]] nabla P_deriv_wrt(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                                const CosineHandle &cosJ, const CosineHandle &cosK) {
  auto firstTerm = derivativeLegendre(cosI.getCos(), A) * cosI.derive_wrt<wrt>() * Legendre(cosJ.getCos(), B) *
                   Legendre(cosK.getCos(), C);
  auto secondTerm = Legendre(cosI.getCos(), A) * derivativeLegendre(cosJ.getCos(), B) * cosJ.derive_wrt<wrt>() *
                    Legendre(cosK.getCos(), C);
  auto thirdTerm = Legendre(cosI.getCos(), A) * Legendre(cosJ.getCos(), B) * derivativeLegendre(cosK.getCos(), C) *
                   cosK.derive_wrt<wrt>();
  return firstTerm + secondTerm + thirdTerm;
}

/**
 *
 * @param A
 * @param B
 * @param C
 * @param cosI
 * @param cosJ
 * @param cosK
 * @return sum of all possible permutations of P(A, B, C, cosI, cosJ, cosK)
 */
[[nodiscard]] double Permutation(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                                 const CosineHandle &cosJ, const CosineHandle &cosK) {
  return P(A, B, C, cosI, cosJ, cosK) + P(A, C, B, cosI, cosJ, cosK) + P(B, A, C, cosI, cosJ, cosK) +
         P(B, C, A, cosI, cosJ, cosK) + P(C, A, B, cosI, cosJ, cosK) + P(C, B, A, cosI, cosJ, cosK);
}

/**
 *
 * @tparam wrt particle with respect to which we are computing the derivative
 * @param A
 * @param B
 * @param C
 * @param cosI
 * @param cosJ
 * @param cosK
 * @return derivative of the permutation summation
 */
template <size_t wrt>
[[nodiscard]] nabla Permutation_deriv_wrt(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                                          const CosineHandle &cosJ, const CosineHandle &cosK) {
  return P_deriv_wrt<wrt>(A, B, C, cosI, cosJ, cosK) + P_deriv_wrt<wrt>(A, C, B, cosI, cosJ, cosK) +
         P_deriv_wrt<wrt>(B, A, C, cosI, cosJ, cosK) + P_deriv_wrt<wrt>(B, C, A, cosI, cosJ, cosK) +
         P_deriv_wrt<wrt>(C, A, B, cosI, cosJ, cosK) + P_deriv_wrt<wrt>(C, B, A, cosI, cosJ, cosK);
}

}  // namespace autopas::utils::ArrayMath::Argon
