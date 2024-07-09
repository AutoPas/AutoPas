/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 13/06/24
 */

#pragma once

#include <cmath>

#include "CosineHandle.h"
#include "DisplacementHandle.h"
#include "Legendre.h"
#include "Parameters.h"

namespace autopas::utils::ArrayMath::Argon {

/**
 *
 * @tparam ID id of the particle with respect to which we are computing the derivative
 * @param a first index of the summation
 * @param b second index of the summation
 * @param c third index of the summation
 * @param A list of values of parameter A
 * @param alpha list of values of parameter alpha
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @return each term of the summation of the repulsive contribution to the force acting on particle with id ID
 */
template <size_t a, size_t b, size_t c, size_t ID>
[[nodiscard]] nabla F_repulsive_abc(const std::array<double, 23> &A, const std::array<double, 23> &alpha,
                                    DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                    DisplacementHandle displacementKI) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto cosineI = CosineHandle(displacementIJ, displacementKI.getInv());
  const auto cosineJ = CosineHandle(displacementIJ.getInv(), displacementJK);
  const auto cosineK = CosineHandle(displacementKI, displacementJK.getInv());

  const auto A_abc = A[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];
  const auto alpha_abc = alpha[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];

  const auto multiplyingFactor{-A_abc * std::exp(-alpha_abc * (IJ + JK + KI))};

  const auto nablaIJ = displacementIJ.derive_wrt<ID>();
  const auto nablaJK = displacementJK.derive_wrt<ID>();
  const auto nablaKI = displacementKI.derive_wrt<ID>();

  const auto firstTerm{ (nablaIJ + nablaJK + nablaJK) * (-alpha_abc)};
  firstTerm = firstTerm * Permutation(a, b, c, cosineI, cosineJ, cosineK);

  const auto secondTerm{Permutation_deriv_wrt<ID>(a, b, c, cosineI, cosineJ, cosineK)};

  return multiplyingFactor * (firstTerm + secondTerm);
}

template <size_t a, size_t b, size_t c>
[[nodiscard]] double U_repulsive_abc(const std::array<double, 23> &A, const std::array<double, 23> &alpha,
                                    DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                    DisplacementHandle displacementKI) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto cosineI = CosineHandle(displacementIJ, displacementKI.getInv());
  const auto cosineJ = CosineHandle(displacementIJ.getInv(), displacementJK);
  const auto cosineK = CosineHandle(displacementKI, displacementJK.getInv());

  const auto A_abc = A[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];
  const auto alpha_abc = alpha[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];

  const auto permutationTerm = Permutation(a, b, c, cosineI, cosineJ, cosineK);

  return A_abc * std::exp(- alpha_abc * (IJ + JK + KI)) * permutationTerm;
}

    /**
 *
 * @tparam ID ID of the particle with respect to whose position we are calculating the derivative
 * @param A list of values of parameter A
 * @param alpha list of values of parameter alpha
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @return Repulsive force acting on the particle with id ID
 */
template <size_t ID>
[[nodiscard]] nabla F_repulsive(const std::array<double, 23> &A, const std::array<double, 23> &alpha,
                                DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                DisplacementHandle displacementKI) {
  const auto F = F_repulsive_abc<0, 0, 0, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 0, 1, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 1, 1, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<1, 1, 1, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 0, 2, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 1, 2, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<1, 1, 2, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 2, 2, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<1, 2, 2, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<2, 2, 2, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 0, 3, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 1, 3, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<1, 1, 3, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 2, 3, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<1, 2, 3, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 3, 3, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 0, 4, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 1, 4, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<1, 1, 4, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 2, 4, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 0, 5, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 1, 5, ID>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 F_repulsive_abc<0, 0, 6, ID>(A, alpha, displacementIJ, displacementJK, displacementKI);
  return F;
}

[[nodiscard]] double U_repulsive(const std::array<double, 23> &A, const std::array<double, 23> &alpha,
                                DisplacementHandle displacementIJ, DisplacementHandle displacementJK,
                                DisplacementHandle displacementKI) {
  const auto U = U_repulsive_abc<0, 0, 0>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 0, 1>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 1, 1>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<1, 1, 1>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 0, 2>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 1, 2>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<1, 1, 2>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 2, 2>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<1, 2, 2>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<2, 2, 2>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 0, 3>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 1, 3>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<1, 1, 3>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 2, 3>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<1, 2, 3>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 3, 3>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 0, 4>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 1, 4>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<1, 1, 4>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 2, 4>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 0, 5>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 1, 5>(A, alpha, displacementIJ, displacementJK, displacementKI) +
                 U_repulsive_abc<0, 0, 6>(A, alpha, displacementIJ, displacementJK, displacementKI);
  return U;
}
}  // namespace autopas::utils::ArrayMath::Argon