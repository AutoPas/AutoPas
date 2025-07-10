//
// Created by sliep on 19.06.2025.
//
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
[[nodiscard]] double RepTerm_abc(const std::array<double, 23> &A, const std::array<double, 23> &alpha,
                                 const double distSquaredIJ, const double distSquaredKI, const double distSquaredJK,
                                 double KIcosIJ, double IJcosJK,
                                 double JKcosKI) {  //        const auto IJ = L2Norm(displacementIJ.getDisplacement());
                                                    //        const auto JK = L2Norm(displacementJK.getDisplacement());
                                                    //        const auto KI = L2Norm(displacementKI.getDisplacement());

  double distIJ = std::sqrt(distSquaredIJ);
  double distJK = std::sqrt(distSquaredJK);
  double distKI = std::sqrt(distSquaredKI);
  /*const auto cosineI = CosineHandle(displacementIJ, displacementKI.getInv());
  const auto cosineJ = CosineHandle(displacementIJ.getInv(), displacementJK);
  const auto cosineK = CosineHandle(displacementKI, displacementJK.getInv());*/

  const auto A_abc = A[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];
  const auto alpha_abc = alpha[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];

  const auto multiplyingFactor{-A_abc * std::exp(-alpha_abc * (distIJ + distJK + distKI))};

  const auto firstTerm{-alpha_abc * (distIJ + distJK + distKI) * Permutation(a, b, c, KIcosIJ, IJcosJK, JKcosKI)};

  //        const auto secondTerm{Permutation_deriv_wrt<ID>(a, b, c, KIcosIJ, IJcosJK, JKcosKI)};

  //        auto res  = multiplyingFactor * (firstTerm + secondTerm);
  //        return multiplyingFactor * (firstTerm + secondTerm);
  return multiplyingFactor * (firstTerm);
}

template <size_t ID>
//    [[nodiscard]] nabla F_RepTermBig(const std::array<double, 23> &A, const std::array<double, 23> &alpha,
[[nodiscard]] double F_RepTermBig(const std::array<double, 23> &A, const std::array<double, 23> &alpha,
                                  const double distSquaredIJ, const double distSquaredKI, const double distSquaredJK) {
  double distIJ = std::sqrt(distSquaredIJ);
  double distJK = std::sqrt(distSquaredJK);
  double distKI = std::sqrt(distSquaredKI);

  double cosineI = (distSquaredIJ + distSquaredKI - distSquaredJK) / (2 * distIJ * distKI);
  double cosineJ = (distSquaredIJ + distSquaredJK - distSquaredKI) / (2 * distIJ * distJK);
  double cosineK = (distSquaredJK + distSquaredKI - distSquaredIJ) / (2 * distJK * distKI);

  const auto F =
      RepTerm_abc<0, 0, 0, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 0, 1, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 1, 1, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<1, 1, 1, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 0, 2, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 1, 2, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<1, 1, 2, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 2, 2, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<1, 2, 2, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<2, 2, 2, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 0, 3, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 1, 3, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<1, 1, 3, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 2, 3, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<1, 2, 3, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 3, 3, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 0, 4, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 1, 4, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<1, 1, 4, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 2, 4, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 0, 5, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 1, 5, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK) +
      RepTerm_abc<0, 0, 6, ID>(A, alpha, distSquaredIJ, distSquaredKI, distSquaredJK, cosineI, cosineJ, cosineK);
  return F;
}

};  // namespace autopas::utils::ArrayMath::Argon

// AUTOPASFUNCTORVALIDATION_REPUSLIVETERMV2_H
