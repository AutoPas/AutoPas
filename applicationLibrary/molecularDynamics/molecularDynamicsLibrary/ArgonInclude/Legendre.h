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
 * @param d order of the Legendre polynomial
 * @note the polynomial is calculated up to order 6, as this is the maximum order needed in the calculation of the Force
 * @return Legendre Polynomial of order d for variable x, i.e. P_d(x)
 */
template <class T>
[[nodiscard]] constexpr T LegendrePol(const T &x, const size_t d) {
  if (d > 6) {
    throw ExceptionHandler::AutoPasException(
        "Argon Simulation: asked to compute Legendre polynomial of order higher than 6.");
  }
  switch (d) {
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
 * @param d order of the Legendre polynomial
 * @return derivative w.r.t. x of the Legendre polynomial of order d for variable x, i.e. dP_d(x)/dx
 */
template <class T>
[[nodiscard]] constexpr T derivativeLegendre(const T &x, size_t d) {
  auto Legendre_order_minus_1 = (d == 0) ? 0 : LegendrePol(x, d - 1);
  auto num{Legendre_order_minus_1 - x * LegendrePol(x, d)};
  return d * num / (1 - Math::pow<2>(x));
}

/**
 *
 * @param A order of the Legendre Polynomial for cosI
 * @param B order of the Legendre Polynomial for cosJ
 * @param C order of the Legendre Polynomial for cosK
 * @param cosI
 * @param cosJ
 * @param cosK
 * @return P_A(cosI) * P_B(cosJ) * P_C(cosK)
 */
[[nodiscard]] double Q(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                       const CosineHandle &cosJ, const CosineHandle &cosK) {
  return LegendrePol(cosI.getCos(), A) * LegendrePol(cosJ.getCos(), B) * LegendrePol(cosK.getCos(), C);
}

/**
 *
 * @tparam ID id of the particle with respect to which we are computing the derivative
 * @param A
 * @param B
 * @param C
 * @param cosI
 * @param cosJ
 * @param cosK
 * @return derivative of Q, i.e. dP_A(cosI)/dcosI * dcosI/dID * P_B(cosJ) * P_C(cosK) +  P_A(cosI) * dP_B(cosJ)/dcosJ * dcosJ/dID * P_C(cosK) +  P_A(cosI) * P_B(cosJ) * dP_C(cosK)/dcosK * dcosK/dID
 */
template <size_t ID>
[[nodiscard]] nabla Q_deriv_wrt(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                                const CosineHandle &cosJ, const CosineHandle &cosK) {
  auto firstTerm = derivativeLegendre(cosI.getCos(), A) * cosI.derive_wrt<ID>() * LegendrePol(cosJ.getCos(), B) *
                   LegendrePol(cosK.getCos(), C);
  auto secondTerm = LegendrePol(cosI.getCos(), A) * derivativeLegendre(cosJ.getCos(), B) * cosJ.derive_wrt<ID>() *
                    LegendrePol(cosK.getCos(), C);
  auto thirdTerm = LegendrePol(cosI.getCos(), A) * LegendrePol(cosJ.getCos(), B) * derivativeLegendre(cosK.getCos(), C) *
                   cosK.derive_wrt<ID>();
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
 * @return sum of the six possible terms that result from the permutation of the interior angles
 */
[[nodiscard]] double Permutation(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                                 const CosineHandle &cosJ, const CosineHandle &cosK) {
  return Q(A, B, C, cosI, cosJ, cosK) + Q(A, C, B, cosI, cosJ, cosK) + Q(B, A, C, cosI, cosJ, cosK) +
         Q(B, C, A, cosI, cosJ, cosK) + Q(C, A, B, cosI, cosJ, cosK) + Q(C, B, A, cosI, cosJ, cosK);
}

/**
 *
 * @tparam ID id of the particle with respect to which we are computing the derivative
 * @param A
 * @param B
 * @param C
 * @param cosI
 * @param cosJ
 * @param cosK
 * @return derivative of the permutation summation
 */
template <size_t ID>
[[nodiscard]] nabla Permutation_deriv_wrt(const size_t A, const size_t B, const size_t C, const CosineHandle &cosI,
                                          const CosineHandle &cosJ, const CosineHandle &cosK) {
  return Q_deriv_wrt<ID>(A, B, C, cosI, cosJ, cosK) + Q_deriv_wrt<ID>(A, C, B, cosI, cosJ, cosK) +
         Q_deriv_wrt<ID>(B, A, C, cosI, cosJ, cosK) + Q_deriv_wrt<ID>(B, C, A, cosI, cosJ, cosK) +
         Q_deriv_wrt<ID>(C, A, B, cosI, cosJ, cosK) + Q_deriv_wrt<ID>(C, B, A, cosI, cosJ, cosK);
}

}  // namespace autopas::utils::ArrayMath::Argon
