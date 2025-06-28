/**
* @file ArgonFunctor.h
* @author I. Angelucci
* @date 18/06/24
*/

#pragma once

#include "DisplacementHandle.h"

namespace autopas::utils::ArrayMath::Argon {

namespace {

[[nodiscard]] size_t factorial(const size_t &n) {
 if (n == 0 || n == 1) return 1;
 return n * factorial(n - 1);
}

template <typename T>
[[nodiscard]] T McLaurianSeries(const size_t &from, const size_t &to, const T &x) {
 T result{0};
 for (size_t n = from; n <= to; n++) {
   result += std::pow(x, n) / factorial(n);
 }
 return result;
}

template <typename T>
[[nodiscard]] T McLaurianSeriesDerivativeTerm(const size_t &from, const size_t &to, const T &r, const double &beta) {
 T result{0};
 auto x = beta * r;
 for (size_t n = from; n <= to; n++) {
   result += std::pow(x, n) / factorial(n) * (beta - n / r);
 }
 return result;
}

}  // namespace

/**
*
* @param beta parameter beta
* @param displacementAB displacement between particle A and B
* @param n power to which r_AB is raised in the corresponding angular function W
* @return damping term D of eq. 9
*/
[[nodiscard]] inline double DampingTerm(const double &beta, const DisplacementHandle &displacementAB, const size_t &n) {
 const auto AB = L2Norm(displacementAB.getDisplacement());
 const auto beta_AB = beta * AB;
 const auto exp = std::exp(-beta_AB);
 const auto series = McLaurianSeries(0, n, beta_AB);
 return 1 - exp * series;
}

/**
*
* @tparam ID id of the particle with respect to which we are computing the derivative
* @param beta parameter beta
* @param displacementAB displacement between particle A and B
* @param n power to which r_AB is raised in the corresponding angular function W
* @return derivative of the damping term D of eq. 9 with respect to the position of particle with id ID
*/
template <size_t ID>
[[nodiscard]] nabla DampingTerm_deriv_wrt(const double &beta, const DisplacementHandle &displacementAB,
                                         const size_t &n) {
 const auto AB = L2Norm(displacementAB.getDisplacement());
 const auto exp = std::exp(-beta * AB);
 const auto series = McLaurianSeriesDerivativeTerm(0, n, AB, beta);
 const auto nablaAB = displacementAB.derive_wrt<ID>();
 return nablaAB * exp * series;
}

}  // namespace autopas::utils::ArrayMath::Argon