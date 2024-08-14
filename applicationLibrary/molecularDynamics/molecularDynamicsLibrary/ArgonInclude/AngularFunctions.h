/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 02/07/24
 */

#pragma once

#include "AngularFunctions.h"

/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 19/06/24
 */

#pragma once

#include "CosineHandle.h"
#include "DisplacementHandle.h"

namespace autopas::utils::ArrayMath::Argon {

namespace {

[[nodiscard]] double cos_2theta(const double &cos_theta) { return 2 * Math::pow<2>(cos_theta) - 1; }

template <size_t ID>
[[nodiscard]] nabla cos_2theta_derive_wrt(const CosineHandle &cosine) {
  const auto cos_deriv = cosine.derive_wrt<ID>();
  return 4 * cosine.getCos() * cos_deriv;
}

[[nodiscard]] double cos_3theta(const double &cos_theta) { return 4 * Math::pow<3>(cos_theta) - 3 * cos_theta; }

template <size_t ID>
[[nodiscard]] nabla cos_3theta_derive_wrt(const CosineHandle &cosine) {
  const auto cos_deriv = cosine.derive_wrt<ID>();
  return (12 * cosine.getCos() * cosine.getCos() - 3) * cos_deriv;
}

[[nodiscard]] double cos_4theta(const double &cos_theta) {
  const auto cos_2theta = 2 * Math::pow<2>(cos_theta) - 1;
  return 2 * Math::pow<2>(cos_2theta) - 1;
}

template <size_t ID>
[[nodiscard]] nabla cos_4theta_derive_wrt(const CosineHandle &cosine) {
  const auto cos_theta = cosine.getCos();
  const auto cos_2theta = 2 * Math::pow<2>(cos_theta) - 1;
  const auto cos_deriv = cosine.derive_wrt<ID>();
  return 16 * cos_2theta * cos_theta * cos_deriv;
}

[[nodiscard]] double sin_theta(const double &cos_theta) { return std::sqrt(1 - Math::pow<2>(cos_theta)); }

template <size_t ID>
[[nodiscard]] nabla sin_theta_derive_wrt(const CosineHandle &cosine) {
  const auto cos_theta = cosine.getCos();
  const auto cos_deriv = cosine.derive_wrt<ID>();
  if (cos_deriv == std::array<double, 3>{{0, 0, 0}}) {
    return cos_deriv;
  }
  return -cos_theta / sin_theta(cos_theta) * cos_deriv;
}

[[nodiscard]] double cos_theta1_minus_theta2(const double &cos_theta1, const double &cos_theta_2) {
  return cos_theta1 * cos_theta_2 + sin_theta(cos_theta1) * sin_theta(cos_theta_2);
}

template <size_t ID>
[[nodiscard]] nabla cos_theta1_minus_theta2_derive_wrt(const CosineHandle &cosine1, const CosineHandle &cosine2) {
  const auto cos_1 = cosine1.getCos();
  const auto cos_1_deriv = cosine1.derive_wrt<ID>();
  const auto cos_2 = cosine2.getCos();
  const auto cos_2_deriv = cosine2.derive_wrt<ID>();
  const auto sin_1 = sin_theta(cosine1.getCos());
  const auto sin_1_deriv = sin_theta_derive_wrt<ID>(cosine1);
  const auto sin_2 = sin_theta(cosine2.getCos());
  const auto sin_2_deriv = sin_theta_derive_wrt<ID>(cosine2);

  return cos_1_deriv * cos_2 + cos_1 * cos_2_deriv + sin_1_deriv * sin_2 + sin_1 * sin_2_deriv;
}

[[nodiscard]] double cos_2_theta1_minus_theta2(const double &cos_theta1, const double &cos_theta_2) {
  const auto cos_theta_1_minus_theta2 = cos_theta1_minus_theta2(cos_theta1, cos_theta_2);
  return cos_2theta(cos_theta_1_minus_theta2);
}

template <size_t ID>
[[nodiscard]] nabla cos_2_theta1_minus_theta2_derive_wrt(const CosineHandle &cosine1, const CosineHandle &cosine2) {
  const auto cos_1 = cosine1.getCos();
  const auto cos_2 = cosine2.getCos();
  const auto cos_theta_1_minus_theta_2 = cos_theta1_minus_theta2(cos_1, cos_2);

  return 4 * cos_theta_1_minus_theta_2 * cos_theta1_minus_theta2_derive_wrt<ID>(cosine1, cosine2);
}

[[nodiscard]] double W_111(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();

  return 3. / (Math::pow<3>(IJ) * Math::pow<3>(JK) * Math::pow<3>(KI)) * (1 + 3 * cosI * cosJ * cosK);
}

[[nodiscard]] double W_112(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();

  const auto cos2K = cos_2theta(cosK);
  const auto cos3K = cos_3theta(cosK);
  const auto cosIminusJ = cos_theta1_minus_theta2(cosI, cosJ);

  return 3. / 16 / (Math::pow<3>(IJ) * Math::pow<4>(JK) * Math::pow<4>(KI)) *
         (9 * cosK - 25 * cos3K + 6 * cosIminusJ * (3 + 5 * cos2K));
}

[[nodiscard]] double W_211(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // no need to call .getInv() on the displacements as we only take their L2Norm
  // permutation of i and k
  return W_112(displacementJK, displacementIJ, displacementKI, cosineK, cosineJ, cosineI);
}

[[nodiscard]] double W_121(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // no need to call .getInv() on the displacements as we only take their L2Norm
  // permutation of j and k
  return W_112(displacementKI, displacementJK, displacementIJ, cosineI, cosineK, cosineJ);
}

[[nodiscard]] double W_122(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();

  const auto cos2I = cos_2theta(cosI);
  const auto cos3I = cos_3theta(cosI);
  const auto cosJminusK = cos_theta1_minus_theta2(cosJ, cosK);
  const auto cos2JminusK = cos_2_theta1_minus_theta2(cosJ, cosK);

  return 15. / 64 / (Math::pow<4>(IJ) * Math::pow<5>(JK) * Math::pow<4>(KI)) *
         (3 * (cosI + 5 * cos3I) + 20 * cosJminusK * (1 - 3 * cos2I) + 70 * cos2JminusK * cosI);
}

[[nodiscard]] double W_212(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // no need to call .getInv() on the displacements as we only take their L2Norm
  // permutation of i and j
  return W_122(displacementIJ, displacementKI, displacementJK, cosineJ, cosineI, cosineK);
}

[[nodiscard]] double W_221(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // no need to call .getInv() on the displacements as we only take their L2Norm
  // permutation of i and k
  return W_122(displacementJK, displacementIJ, displacementKI, cosineK, cosineJ, cosineI);
}

[[nodiscard]] double W_222(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();

  const auto cos2I = cos_2theta(cosI);
  const auto cos2J = cos_2theta(cosJ);
  const auto cos2K = cos_2theta(cosK);

  const auto cos2IminusJ = cos_2_theta1_minus_theta2(cosI, cosJ);
  const auto cos2JminusK = cos_2_theta1_minus_theta2(cosJ, cosK);
  const auto cos2KminusI = cos_2_theta1_minus_theta2(cosK, cosI);

  return 15. / 128 / (Math::pow<5>(IJ) * Math::pow<5>(JK) * Math::pow<5>(KI)) *
         (-27 + 220 * cosI * cosJ * cosK + 490 * cos2I * cos2J * cos2K +
          175 * (cos2IminusJ + cos2JminusK + cos2KminusI));
}

[[nodiscard]] double W_113(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();

  const auto cos2K = cos_2theta(cosK);
  const auto cos3K = cos_3theta(cosK);
  const auto cos4K = cos_4theta(cosK);
  const auto cosIminusJ = cos_theta1_minus_theta2(cosI, cosJ);

  return 5. / 32 / (Math::pow<3>(IJ) * Math::pow<5>(JK) * Math::pow<5>(KI)) *
         (9 + 8 * cos2K - 49 * cos4K + 6 * cosIminusJ * (9 * cosK + 7 * cos3K));
}

[[nodiscard]] double W_311(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // no need to call .getInv() on the displacements as we only take their L2Norm
  // permutation of i and k
  return W_113(displacementJK, displacementIJ, displacementKI, cosineK, cosineJ, cosineI);
}

[[nodiscard]] double W_131(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // no need to call .getInv() on the displacements as we only take their L2Norm
  // permutation of j and k
  return W_113(displacementKI, displacementJK, displacementIJ, cosineI, cosineK, cosineJ);
}

template <size_t ID>
[[nodiscard]] nabla W_111_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto nablaIJ = displacementIJ.derive_wrt<ID>();
  const auto nablaJK = displacementJK.derive_wrt<ID>();
  const auto nablaKI = displacementKI.derive_wrt<ID>();

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();

  const auto nablaCosI = cosineI.derive_wrt<ID>();
  const auto nablaCosJ = cosineJ.derive_wrt<ID>();
  const auto nablaCosK = cosineK.derive_wrt<ID>();

  const auto nablaDisplacementTerm = -1. * (nablaIJ / IJ + nablaJK / JK + nablaKI / KI);
  const auto cosineTerm = 1 + 3 * cosI * cosJ * cosK;
  const auto nablaCosineTerm = nablaCosI * cosJ * cosK + cosI * nablaCosJ * cosK + cosI * cosJ * nablaCosK;

  return 9. / (Math::pow<3>(IJ) * Math::pow<3>(JK) * Math::pow<3>(KI)) *
         (nablaDisplacementTerm * cosineTerm + nablaCosineTerm);
}

template <size_t ID>
[[nodiscard]] nabla W_112_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto nablaIJ = displacementIJ.derive_wrt<ID>();
  const auto nablaJK = displacementJK.derive_wrt<ID>();
  const auto nablaKI = displacementKI.derive_wrt<ID>();

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();
  const auto cos2K = cos_2theta(cosK);
  const auto cos3K = cos_3theta(cosK);
  const auto cosIminusJ = cos_theta1_minus_theta2(cosI, cosJ);

  const auto nablaCosK = cosineK.derive_wrt<ID>();
  const auto nablaCos2K = cos_2theta_derive_wrt<ID>(cosineK);
  const auto nablaCos3K = cos_3theta_derive_wrt<ID>(cosineK);
  const auto nablaCosIminusJ = cos_theta1_minus_theta2_derive_wrt<ID>(cosineI, cosineJ);

  const auto nablaDisplacementTerm = -3. * nablaIJ / IJ - 4. * nablaJK / JK - 4. * nablaKI / KI;
  const auto cosineTerm = 9. * cosK - 25. * cos3K + 6. * cosIminusJ * (3 + 5. * cos2K);
  const auto nablaCosineTerm =
      9. * nablaCosK - 25. * nablaCos3K + 6. * nablaCosIminusJ * (3 + 5. * cos2K) + 30. * cosIminusJ * nablaCos2K;

  return 3. / 16 / (Math::pow<3>(IJ) * Math::pow<4>(JK) * Math::pow<4>(KI)) *
         (nablaDisplacementTerm * cosineTerm + nablaCosineTerm);
}

template <size_t ID>
[[nodiscard]] nabla W_211_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // permutation of i and k
  return W_112_deriv_wrt<ID>(displacementJK.getInv(), displacementIJ.getInv(), displacementKI.getInv(), cosineK,
                             cosineJ, cosineI);
}

template <size_t ID>
[[nodiscard]] nabla W_121_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // permutation of j and k
  return W_112_deriv_wrt<ID>(displacementKI.getInv(), displacementJK.getInv(), displacementIJ.getInv(), cosineI,
                             cosineK, cosineJ);
}

template <size_t ID>
[[nodiscard]] nabla W_122_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto nablaIJ = displacementIJ.derive_wrt<ID>();
  const auto nablaJK = displacementJK.derive_wrt<ID>();
  const auto nablaKI = displacementKI.derive_wrt<ID>();

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();
  const auto cos2I = cos_2theta(cosI);
  const auto cos3I = cos_3theta(cosI);
  const auto cosJminusK = cos_theta1_minus_theta2(cosJ, cosK);
  const auto cos2JminusK = cos_2_theta1_minus_theta2(cosJ, cosK);

  const auto nablaCosI = cosineI.derive_wrt<ID>();
  const auto nablaCos2I = cos_2theta_derive_wrt<ID>(cosineI);
  const auto nablaCos3I = cos_3theta_derive_wrt<ID>(cosineI);
  const auto nablaCosJminusK = cos_theta1_minus_theta2_derive_wrt<ID>(cosineJ, cosineK);
  const auto nablaCos2JminusK = cos_2_theta1_minus_theta2_derive_wrt<ID>(cosineJ, cosineK);

  const auto nablaDisplacementTerm = -4. * nablaIJ / IJ - 5. * nablaJK / JK - 4. * nablaKI / KI;
  const auto cosineTerm = 3. * cosI + 15. * cos3I + 20. * cosJminusK * (1 - 3. * cos2I) + 70. * cos2JminusK * cosI;
  const auto nablaCosineTerm = 3. * nablaCosI + 15. * nablaCos3I + 20. * nablaCosJminusK * (1 - 3. * cos2I) -
                               60. * cosJminusK * nablaCos2I + 70. * nablaCos2JminusK * cosI +
                               70. * cos2JminusK * nablaCosI;

  return 15. / 64 / (Math::pow<4>(IJ) * Math::pow<5>(JK) * Math::pow<4>(KI)) *
         (nablaDisplacementTerm * cosineTerm + nablaCosineTerm);
}

template <size_t ID>
[[nodiscard]] nabla W_212_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // permutation of i and j
  return W_122_deriv_wrt<ID>(displacementIJ.getInv(), displacementKI.getInv(), displacementJK.getInv(), cosineJ,
                             cosineI, cosineK);
}

template <size_t ID>
[[nodiscard]] nabla W_221_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // permutation of i and k
  return W_122_deriv_wrt<ID>(displacementJK.getInv(), displacementIJ.getInv(), displacementKI.getInv(), cosineK,
                             cosineJ, cosineI);
}

template <size_t ID>
[[nodiscard]] nabla W_222_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto nablaIJ = displacementIJ.derive_wrt<ID>();
  const auto nablaJK = displacementJK.derive_wrt<ID>();
  const auto nablaKI = displacementKI.derive_wrt<ID>();

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();
  const auto cos2I = cos_2theta(cosI);
  const auto cos2J = cos_2theta(cosJ);
  const auto cos2K = cos_2theta(cosK);
  const auto cos2IminusJ = cos_2_theta1_minus_theta2(cosI, cosJ);
  const auto cos2JminusK = cos_2_theta1_minus_theta2(cosJ, cosK);
  const auto cos2KminusI = cos_2_theta1_minus_theta2(cosK, cosI);

  const auto nablaCosI = cosineI.derive_wrt<ID>();
  const auto nablaCosJ = cosineJ.derive_wrt<ID>();
  const auto nablaCosK = cosineK.derive_wrt<ID>();
  const auto nablaCos2I = cos_2theta_derive_wrt<ID>(cosineI);
  const auto nablaCos2J = cos_2theta_derive_wrt<ID>(cosineJ);
  const auto nablaCos2K = cos_2theta_derive_wrt<ID>(cosineK);
  const auto nablaCos2IminusJ = cos_2_theta1_minus_theta2_derive_wrt<ID>(cosineI, cosineJ);
  const auto nablaCos2JminusK = cos_2_theta1_minus_theta2_derive_wrt<ID>(cosineJ, cosineK);
  const auto nablaCos2KminusI = cos_2_theta1_minus_theta2_derive_wrt<ID>(cosineK, cosineI);

  const auto nablaDisplacementTerm = -5. * (nablaIJ / IJ + nablaJK / JK + nablaKI / KI);
  const auto cosineTerm = -27. + 220. * cosI * cosJ * cosK + 490. * cos2I * cos2J * cos2K +
                          175. * (cos2IminusJ + cos2JminusK + cos2KminusI);
  const auto nablaCosineTerm1 = 220. * (nablaCosI * cosJ * cosK + cosI * nablaCosJ * cosK + cosI * cosJ * nablaCosK);
  const auto nablaCosineTerm2 =
      490. * (nablaCos2I * cos2J * cos2K + cos2I * nablaCos2J * cos2K + cos2I * cos2J * nablaCos2K);
  const auto nablaCosineTerm3 = 175. * (nablaCos2IminusJ + nablaCos2JminusK + nablaCos2KminusI);

  return 15. / 128 / (Math::pow<5>(IJ) * Math::pow<5>(JK) * Math::pow<5>(KI)) *
         (nablaDisplacementTerm * cosineTerm + nablaCosineTerm1 + nablaCosineTerm2 + nablaCosineTerm3);
}

template <size_t ID>
[[nodiscard]] nabla W_113_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  const auto IJ = L2Norm(displacementIJ.getDisplacement());
  const auto JK = L2Norm(displacementJK.getDisplacement());
  const auto KI = L2Norm(displacementKI.getDisplacement());

  const auto nablaIJ = displacementIJ.derive_wrt<ID>();
  const auto nablaJK = displacementJK.derive_wrt<ID>();
  const auto nablaKI = displacementKI.derive_wrt<ID>();

  const auto cosI = cosineI.getCos();
  const auto cosJ = cosineJ.getCos();
  const auto cosK = cosineK.getCos();
  const auto cos2K = cos_2theta(cosK);
  const auto cos3K = cos_3theta(cosK);
  const auto cos4K = cos_4theta(cosK);
  const auto cosIminusJ = cos_theta1_minus_theta2(cosI, cosJ);

  const auto nablaCosK = cosineK.derive_wrt<ID>();
  const auto nablaCos2K = cos_2theta_derive_wrt<ID>(cosineK);
  const auto nablaCos3K = cos_3theta_derive_wrt<ID>(cosineK);
  const auto nablaCos4K = cos_4theta_derive_wrt<ID>(cosineK);
  const auto nablaCosIminusJ = cos_theta1_minus_theta2_derive_wrt<ID>(cosineI, cosineJ);

  const auto nablaDisplacementTerm = -3. * nablaIJ / IJ - 5. * nablaJK / JK - 5. * nablaKI / KI;
  const auto cosineTerm = 9. + 8. * cos2K - 49. * cos4K + 6. * cosIminusJ * (9. * cosK + 7. * cos3K);
  const auto nablaCosineTerm = 8. * nablaCos2K - 49. * nablaCos4K + 6. * nablaCosIminusJ * (9. * cosK + 7. * cos3K) +
                               6. * cosIminusJ * (9. * nablaCosK + 7. * nablaCos3K);

  return 5. / 32 / (Math::pow<3>(IJ) * Math::pow<5>(JK) * Math::pow<5>(KI)) *
         (nablaDisplacementTerm * cosineTerm + nablaCosineTerm);
}

template <size_t ID>
[[nodiscard]] nabla W_131_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // permuation of j and k
  return W_113_deriv_wrt<ID>(displacementKI.getInv(), displacementJK.getInv(), displacementIJ.getInv(), cosineI,
                             cosineK, cosineJ);
}

template <size_t ID>
[[nodiscard]] nabla W_311_deriv_wrt(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                    const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                    const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  // permutation of i and k
  return W_113_deriv_wrt<ID>(displacementJK.getInv(), displacementIJ.getInv(), displacementKI.getInv(), cosineK,
                             cosineJ, cosineI);
}

}  // namespace

/**
 *
 * @tparam a first l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam b second l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam c third l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @param cosineI angle between the displacements IJ and IK
 * @param cosineJ angle between the displacements JK and JI
 * @param cosineK angle between the displacements KI and KJ
 * @return angular function derived by Bell describing the interaction between dipole/quadrupole/octupole
 */
template <size_t a, size_t b, size_t c>
[[nodiscard]] double AngularTerm(const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                 const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                 const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  if (a == 1 && b == 1 && c == 1) {
    return W_111(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 1 && c == 2) {
    return W_112(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 2 && c == 1) {
    return W_121(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 1 && c == 1) {
    return W_211(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 2 && c == 2) {
    return W_122(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 1 && c == 2) {
    return W_212(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 2 && c == 1) {
    return W_221(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 2 && c == 2) {
    return W_222(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 1 && c == 3) {
    return W_113(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 3 && c == 1) {
    return W_131(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 3 && b == 1 && c == 1) {
    return W_311(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
}

/**
 *
 * @tparam a first l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam b second l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam c third l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam ID ID of the particle with respect to whose position we are calculating the derivative
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @param cosineI angle between the displacements IJ and IK
 * @param cosineJ angle between the displacements JK and JI
 * @param cosineK angle between the displacements KI and KJ
 * @return derivative of angular function derived by Bell describing the interaction between dipole/quadrupole/octupole
 */
template <size_t a, size_t b, size_t c, size_t ID>
[[nodiscard]] nabla AngularTerm_derive_wrt(const DisplacementHandle &displacementIJ,
                                           const DisplacementHandle &displacementJK,
                                           const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                           const CosineHandle &cosineJ, const CosineHandle &cosineK) {
  if (a == 1 && b == 1 && c == 1) {
    return W_111_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 1 && c == 2) {
    return W_112_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 2 && c == 1) {
    return W_121_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 1 && c == 1) {
    return W_211_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 2 && c == 2) {
    return W_122_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 1 && c == 2) {
    return W_212_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 2 && c == 1) {
    return W_221_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 2 && b == 2 && c == 2) {
    return W_222_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 1 && c == 3) {
    return W_113_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 1 && b == 3 && c == 1) {
    return W_131_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
  if (a == 3 && b == 1 && c == 1) {
    return W_311_deriv_wrt<ID>(displacementIJ, displacementJK, displacementKI, cosineI, cosineJ, cosineK);
  }
}

}  //  namespace autopas::utils::ArrayMath::Argon
