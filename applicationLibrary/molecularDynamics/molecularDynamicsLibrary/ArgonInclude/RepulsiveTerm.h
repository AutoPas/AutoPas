/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 13/06/24
 */

#pragma once

#include "CosineHandle.h"
#include "DisplacementHandle.h"
#include "Legendre.h"
#include "Parameters.h"
#include <cmath>

namespace autopas::utils::ArrayMath::Argon {

  /**
   *
   * @tparam particle particle on which the repulsive force is being calculated on
   * @param a first index of the summation
   * @param b second index of the summation
   * @param c third index of the summation
   * @param A parametr A
   * @param alpha parameter alpha
   * @param displacementIJ
   * @param displacementJK
   * @param displacementKI
   * @return single term in the summation of the repulsive force
   */
  template<size_t particle>
  [[nodiscard]] nabla F_repulsive_abc(const size_t& a, const size_t& b, const size_t& c, const std::array<double, 23>& A, const std::array<double, 23>& alpha,
                          DisplacementHandle displacementIJ, DisplacementHandle displacementJK, DisplacementHandle displacementKI) {
    const auto IJ = displacementIJ.getDisplacement();
    const auto JK = displacementJK.getDisplacement();
    const auto KI = displacementKI.getDisplacement();

    const auto cosineI = CosineHandle(displacementIJ, displacementKI.getInv());
    const auto cosineJ = CosineHandle(displacementIJ.getInv(), displacementJK);
    const auto cosineK = CosineHandle(displacementKI, displacementJK.getInv());

    const auto A_abc= A[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];
    const auto alpha_abc= alpha[mdLib::Argon::index<mdLib::Argon::param::A>(a, b, c)];

    const auto multiplyingFactor{-A_abc * std::exp(-alpha_abc * (L2Norm(IJ) + L2Norm(JK) + L2Norm(KI)))};

    const auto firstTerm{-alpha_abc * (displacementIJ.derive_wrt<particle>() + displacementJK.derive_wrt<particle>() + displacementKI.derive_wrt<particle>())};
    firstTerm *= Permutation(a, b, c, cosineI, cosineJ, cosineK);

    const auto secondTerm{Permutation_deriv_wrt<particle>(a, b, c, cosineI, cosineJ, cosineK)};

    return multiplyingFactor * (firstTerm + secondTerm);
  }

  /**
   *
   * @tparam index index on which the repulsive force is being calculated on
   * @param A parameter A
   * @param alpha parameter alpha
   * @param displacementIJ
   * @param displacementJK
   * @param displacementKI
   * @return Repulsive force acting on the index
   */
  template<size_t particle>
  [[nodiscard]] nabla F_repulsive(const std::array<double, 23>& A, const std::array<double, 23>& alpha, DisplacementHandle displacementIJ, DisplacementHandle displacementJK, DisplacementHandle displacementKI) {
    const auto F = F_repulsive_abc<particle>(0, 0, 0, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 0, 1, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 1, 1, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(1, 1, 1, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 0, 2, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 1, 2, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(1, 1, 2, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 2, 2, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(1, 2, 2, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(2, 2, 2, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 0, 3, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 1, 3, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(1, 1, 3, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 2, 3, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(1, 2, 3, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 3, 3, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 0, 4, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 1, 4, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(1, 1, 4, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 2, 4, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 0, 5, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 1, 5, A, alpha, displacementIJ, displacementJK, displacementKI) +
                   F_repulsive_abc<particle>(0, 0, 6, A, alpha, displacementIJ, displacementJK, displacementKI);
    return F;
  }
}   // namespace autopas::utils::ArrayMath::Argon