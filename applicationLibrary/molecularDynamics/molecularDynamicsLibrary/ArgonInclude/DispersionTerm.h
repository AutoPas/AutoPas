/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 02/07/24
 */

#pragma once

#include <cmath>

#include "AngularFunctions.h"
#include "DampingTerm.h"
#include "Parameters.h"

namespace autopas::utils::ArrayMath::Argon {

/**
 *
 * @tparam ID id of the particle with respect to which we are computing the derivative
 * @param a first index of the summation
 * @param b second index of the summation
 * @param c third index of the summation
 * @param Z list of values of parameter Z
 * @param beta list of values of parameter beta
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @return each term of the summation of the dispersion contribution to the force acting on particle with id ID
 */
template <size_t a, size_t b, size_t c, size_t ID>
[[nodiscard]] nabla F_dispersive_abc(const std::array<double, 5> &Z, const std::array<double, 5> &beta,
                                     DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                     DisplacementHandle displacementKI) {
  const auto Z_abc = Z[mdLib::Argon::index<mdLib::Argon::param::Z>(a, b, c)];
  const auto beta_abc = beta[mdLib::Argon::index<mdLib::Argon::param::beta>(a, b, c)];

  const auto cosineI = CosineHandle(displacementIJ, displacementKI.getInv());
  const auto cosineJ = CosineHandle(displacementIJ.getInv(), displacementJK);
  const auto cosineK = CosineHandle(displacementKI, displacementJK.getInv());

  const auto n1 = a + b + 1;
  const auto n2 = b + c + 1;
  const auto n3 = c + a + 1;

  const auto D1 = DampingTerm(beta_abc, displacementIJ, n1);
  const auto D2 = DampingTerm(beta_abc, displacementJK, n2);
  const auto D3 = DampingTerm(beta_abc, displacementKI, n3);

  const auto nabla_D1 = DampingTerm_deriv_wrt<ID>(beta_abc, displacementIJ, n1);
  const auto nabla_D2 = DampingTerm_deriv_wrt<ID>(beta_abc, displacementJK, n2);
  const auto nabla_D3 = DampingTerm_deriv_wrt<ID>(beta_abc, displacementKI, n3);

  const auto W = AngularTerm<a, b, c>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);

  const auto nabla_W =
      AngularTerm_derive_wrt<a, b, c, ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);

  return Z_abc * (W * (nabla_D1 * D2 * D3 + D1 * nabla_D2 * D3 + D1 * D2 * nabla_D3) + D1 * D2 * D3 * nabla_W);
}

template <size_t a, size_t b, size_t c>
[[nodiscard]] double U_dispersive_abc(const std::array<double, 5> &Z, const std::array<double, 5> &beta,
                                      DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                      DisplacementHandle displacementKI) {
  const auto Z_abc = Z[mdLib::Argon::index<mdLib::Argon::param::Z>(a, b, c)];
  const auto beta_abc = beta[mdLib::Argon::index<mdLib::Argon::param::beta>(a, b, c)];

  const auto cosineI = CosineHandle(displacementIJ, displacementKI.getInv());
  const auto cosineJ = CosineHandle(displacementIJ.getInv(), displacementJK);
  const auto cosineK = CosineHandle(displacementKI, displacementJK.getInv());

  const auto n1 = a + b + 1;
  const auto n2 = b + c + 1;
  const auto n3 = c + a + 1;

  const auto D1 = DampingTerm(beta_abc, displacementIJ, n1);
  const auto D2 = DampingTerm(beta_abc, displacementJK, n2);
  const auto D3 = DampingTerm(beta_abc, displacementKI, n3);

  const auto W = AngularTerm<a, b, c>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);

  const auto U = D1 * D2 * D3 * W * Z_abc;
  return U;
}

/**
 *
 * @tparam ID ID of the particle with respect to whose position we are calculating the derivative
 * @param Z list of values of parameter Z
 * @param beta list of values of parameter beta
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @return Repulsive force acting on the particle with id ID
 */
template <size_t ID>
[[nodiscard]] nabla F_dispersive(const std::array<double, 5> &Z, const std::array<double, 5> &beta,
                                 DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                 DisplacementHandle displacementKI) {
  const auto F = F_dispersive_abc<1, 1, 1, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<1, 1, 2, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<1, 2, 1, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<2, 1, 1, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<1, 2, 2, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<2, 1, 2, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<2, 2, 1, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<2, 2, 2, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<1, 1, 3, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<1, 3, 1, ID>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 F_dispersive_abc<3, 1, 1, ID>(Z, beta, displacementIJ, displacementJK, displacementKI);
  return F;
}

[[nodiscard]] double U_dispersive(const std::array<double, 5> &Z, const std::array<double, 5> &beta,
                                  DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                  DisplacementHandle displacementKI) {
  const auto U = U_dispersive_abc<1, 1, 1>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<1, 1, 2>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<1, 2, 1>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<2, 1, 1>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<1, 2, 2>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<2, 1, 2>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<2, 2, 1>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<2, 2, 2>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<1, 1, 3>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<1, 3, 1>(Z, beta, displacementIJ, displacementJK, displacementKI) +
                 U_dispersive_abc<3, 1, 1>(Z, beta, displacementIJ, displacementJK, displacementKI);
  return U;
}

}  // namespace autopas::utils::ArrayMath::Argon
