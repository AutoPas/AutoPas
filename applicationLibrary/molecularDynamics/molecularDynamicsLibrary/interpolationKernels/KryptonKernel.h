/**
 * @file KryptonKernel.h
 *
 * @date 29.04.2025
 * @author Luis Gall
 */

#pragma once

#include "Kernel.h"
#include "autopas/utils/ArrayMath.h"

namespace mdLib {

class KryptonKernel : public Kernel<KryptonKernel> {
 public:
  explicit KryptonKernel()
      : Kernel<KryptonKernel>(){

        };

  double calculatePairDerivative(double dr) final {
    const double dr2 = dr * dr;
    const double distInv = 1. / dr;
    const double distInv2 = distInv * distInv;
    const double distInv6 = distInv2 * distInv2 * distInv2;
    const double distNeg8 = distInv6 * distInv2;
    const double distNeg10 = distNeg8 * distInv2;
    const double distNeg12 = distNeg10 * distInv2;
    const double distNeg14 = distNeg12 * distInv2;
    const double distNeg16 = distNeg14 * distInv2;

    const double expAlphaTerm = std::exp(_alpha1 * dr + _alpha2 * dr2 + _alphaInv1 * distInv);
    const double alphaTerm = _alpha1 + 2 * _alpha2 * dr - _alphaInv1 * distInv2;

    const double firstTerm = (-_constA * alphaTerm * expAlphaTerm) * distInv;

    const double bdist = _constb * dr;
    const double bdist2 = bdist * bdist;
    const double bdist3 = bdist2 * bdist;
    const double bdist4 = bdist3 * bdist;
    const double bdist5 = bdist4 * bdist;
    const double bdist6 = bdist5 * bdist;
    const double bdist7 = bdist6 * bdist;
    const double bdist8 = bdist7 * bdist;
    const double bdist9 = bdist8 * bdist;
    const double bdist10 = bdist9 * bdist;
    const double bdist11 = bdist10 * bdist;
    const double bdist12 = bdist11 * bdist;
    const double bdist13 = bdist12 * bdist;
    const double bdist14 = bdist13 * bdist;
    const double bdist15 = bdist14 * bdist;
    const double bdist16 = bdist15 * bdist;

    const double ksum0 = 1.0;
    const double ksum1 = bdist;
    const double ksum2 = bdist2 * _invFactorials[2];
    const double ksum3 = bdist3 * _invFactorials[3];
    const double ksum4 = bdist4 * _invFactorials[4];
    const double ksum5 = bdist5 * _invFactorials[5];
    const double ksum6 = bdist6 * _invFactorials[6];
    const double ksum7 = bdist7 * _invFactorials[7];
    const double ksum8 = bdist8 * _invFactorials[8];
    const double ksum9 = bdist9 * _invFactorials[9];
    const double ksum10 = bdist10 * _invFactorials[10];
    const double ksum11 = bdist11 * _invFactorials[11];
    const double ksum12 = bdist12 * _invFactorials[12];
    const double ksum13 = bdist13 * _invFactorials[13];
    const double ksum14 = bdist14 * _invFactorials[14];
    const double ksum15 = bdist15 * _invFactorials[15];
    const double ksum16 = bdist16 * _invFactorials[16];

    const double ksumacc6 = ksum0 + ksum1 + ksum2 + ksum3 + ksum4 + ksum5 + ksum6;
    const double ksumacc8 = ksumacc6 + ksum7 + ksum8;
    const double ksumacc10 = ksumacc8 + ksum9 + ksum10;
    const double ksumacc12 = ksumacc10 + ksum11 + ksum12;
    const double ksumacc14 = ksumacc12 + ksum13 + ksum14;
    const double ksumacc16 = ksumacc14 + ksum15 + ksum16;

    const double expbr = std::exp(-bdist);

    const double term6 = _constC6 * distInv6 * (-6.0 + expbr * (6.0 * ksumacc6 + bdist * ksum6));
    const double term8 = _constC8 * distNeg8 * (-8.0 + expbr * (8.0 * ksumacc8 + bdist * ksum8));
    const double term10 = _constC10 * distNeg10 * (-10.0 + expbr * (10.0 * ksumacc10 + bdist * ksum10));
    const double term12 = _constC12 * distNeg12 * (-12.0 + expbr * (12.0 * ksumacc12 + bdist * ksum12));
    const double term14 = _constC14 * distNeg14 * (-14.0 + expbr * (14.0 * ksumacc14 + bdist * ksum14));
    const double term16 = _constC16 * distNeg16 * (-16.0 + expbr * (16.0 * ksumacc16 + bdist * ksum16));

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16) * distInv2;

    return firstTerm + secondTerm;
  }

  double calculatePair(double dr) final {
    const double dr2 = dr * dr;
    const double distInv = 1. / dr;
    const double distInv2 = distInv * distInv;
    const double distInv6 = distInv2 * distInv2 * distInv2;
    const double distNeg8 = distInv6 * distInv2;
    const double distNeg10 = distNeg8 * distInv2;
    const double distNeg12 = distNeg10 * distInv2;
    const double distNeg14 = distNeg12 * distInv2;
    const double distNeg16 = distNeg14 * distInv2;

    const double firstTerm = _constA * std::exp(_alpha1*dr + _alpha2*dr2 + _alphaInv1*distInv);

    const double bdist = _constb * dr;
    const double bdist2 = bdist * bdist;
    const double bdist3 = bdist2 * bdist;
    const double bdist4 = bdist3 * bdist;
    const double bdist5 = bdist4 * bdist;
    const double bdist6 = bdist5 * bdist;
    const double bdist7 = bdist6 * bdist;
    const double bdist8 = bdist7 * bdist;
    const double bdist9 = bdist8 * bdist;
    const double bdist10 = bdist9 * bdist;
    const double bdist11 = bdist10 * bdist;
    const double bdist12 = bdist11 * bdist;
    const double bdist13 = bdist12 * bdist;
    const double bdist14 = bdist13 * bdist;
    const double bdist15 = bdist14 * bdist;
    const double bdist16 = bdist15 * bdist;

    const double ksum0 = 1.0;
    const double ksum1 = bdist;
    const double ksum2 = bdist2 * _invFactorials[2];
    const double ksum3 = bdist3 * _invFactorials[3];
    const double ksum4 = bdist4 * _invFactorials[4];
    const double ksum5 = bdist5 * _invFactorials[5];
    const double ksum6 = bdist6 * _invFactorials[6];
    const double ksum7 = bdist7 * _invFactorials[7];
    const double ksum8 = bdist8 * _invFactorials[8];
    const double ksum9 = bdist9 * _invFactorials[9];
    const double ksum10 = bdist10 * _invFactorials[10];
    const double ksum11 = bdist11 * _invFactorials[11];
    const double ksum12 = bdist12 * _invFactorials[12];
    const double ksum13 = bdist13 * _invFactorials[13];
    const double ksum14 = bdist14 * _invFactorials[14];
    const double ksum15 = bdist15 * _invFactorials[15];
    const double ksum16 = bdist16 * _invFactorials[16];

    const double ksumacc6 = ksum0 + ksum1 + ksum2 + ksum3 + ksum4 + ksum5 + ksum6;
    const double ksumacc8 = ksumacc6 + ksum7 + ksum8;
    const double ksumacc10 = ksumacc8 + ksum9 + ksum10;
    const double ksumacc12 = ksumacc10 + ksum11 + ksum12;
    const double ksumacc14 = ksumacc12 + ksum13 + ksum14;
    const double ksumacc16 = ksumacc14 + ksum15 + ksum16;

    const double expbr = std::exp(-bdist);

    const double term6 = _constC6 * distInv6 * (1.-expbr * ksumacc6);
    const double term8 = _constC8 * distNeg8 * (1.-expbr * ksumacc8);
    const double term10 = _constC10 * distNeg10 * (1.-expbr * ksumacc10);
    const double term12 = _constC12 * distNeg12 * (1.-expbr * ksumacc12);
    const double term14 = _constC14 * distNeg14 * (1.-expbr * ksumacc14);
    const double term16 = _constC16 * distNeg16 * (1.-expbr * ksumacc16);

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16);

    return firstTerm - secondTerm;
  }

  double calculateTriplet(double dr1, double dr2, double dr3) {

    const double firstCos = dr1*dr1 + dr2*dr2 - dr3*dr3;
    const double secondCos = dr1*dr1 - dr2*dr2 + dr3*dr3;
    const double thirdCos = -dr1*dr1 + dr2*dr2 + dr3*dr3;

    const double firstNom = firstCos * secondCos * thirdCos;

    const double firstDenom = dr1 * dr2 * dr3;
    const double firstDenom2 = firstDenom * firstDenom;

    const double firstPart = 1 + (3.*firstNom) / (8.*firstDenom2);

    const double secondDenom3 = firstDenom2 * firstDenom;
    const double exponent = -_alpha * (dr1 + dr2 + dr3);
    const double ePart = std::exp(exponent);

    double sum = 0.;
    for (int i = 0; i <= 5; ++i) {
      sum += _A.at(i) * std::pow(firstDenom, 2.*i/3.);
    }

    const double secondPart = _constC / secondDenom3 + ePart * sum;

    return firstPart * secondPart;
  }

  std::array<double, 3> calculateTripletDerivative(double drIJ, double drJK, double drKI) {
      // see Markus Branch for reference: https://github.com/AutoPas/AutoPas/blob/feat/3xa/noble-gas-functors/applicationLibrary/molecularDynamics/molecularDynamicsLibrary/KryptonExtendedATMFunctor.h
      // dr1 = dIJ
      // dr2 = dJK
      // dr3 = dKI

      const double drIJSquare = drIJ*drIJ;
      const double drJKSquare = drJK*drJK;
      const double drKISquare = drKI*drKI;

      const double numeratorKI = drIJSquare + drJKSquare - drKISquare;
      const double numeratorJK = drIJSquare + drKISquare - drJKSquare;
      const double numeratorIJ = drJKSquare + drKISquare - drIJSquare;
    
      const double numerator = numeratorKI * numeratorJK * numeratorIJ;

      const double allDrSquare = drIJSquare * drJKSquare * drKISquare;
      const double allDr = drIJ * drJK * drKI;
      const double allDrTripled = allDrSquare * allDr;

      const double allDrTripledGradientIJ = 3. / (allDrTripled * drIJSquare);
      const double allDrTripledGradientKI = -3. / (allDrTripled * drKISquare);

      const double cosines = (3. / 8.) * numerator / allDrSquare;
      const double cosinesGradientIJ = (3./4.) * ((numerator / drIJSquare - numeratorKI * numeratorIJ - numeratorJK * numeratorIJ + numeratorJK * numeratorKI) / allDrSquare);
      const double cosinesGradientKI = (3./4.) * ((-numerator / drKISquare + numeratorKI * numeratorIJ - numeratorJK * numeratorIJ + numeratorJK * numeratorKI) / allDrSquare);
  
      const double fullAtmGradientIJ = _constC * ((1.+cosines) * allDrTripledGradientIJ + cosinesGradientIJ / allDrTripled);
      const double fullAtmGradientKI = _constC * ((1.+cosines) * allDrTripledGradientKI + cosinesGradientKI / allDrTripled);

      const double expTerm = std::exp(-_alpha * (drIJ + drJK + drKI));

      std::array<double, 6> sumFactors {};
      double sum = 0.;
      for (size_t n = 0; n < sumFactors.size(); ++n) {
        sumFactors.at(n) = _A.at(n) * std::pow(drIJ*drJK*drKI, 2.*n / 3.);
        sum += sumFactors.at(n);
      }


      double ijSum = 0.;
      for (size_t n = 0; n < sumFactors.size(); ++n) {
        ijSum += sumFactors.at(n) * (2.*n / (3.*drIJ) - _alpha);
      }

      double kiSum = 0.;
      for (size_t n = 0; n < sumFactors.size(); ++n) {
        kiSum += sumFactors.at(n) * (2.*n / (3.*drKI) - _alpha);
      }

      const double fullExpGradientIJ = expTerm * (-(1.+cosines) * ijSum / drIJ + cosinesGradientIJ * sum);
      const double fullExpGradientKI = expTerm * ((1.+cosines) * kiSum / drKI + cosinesGradientKI * sum);

      const double gradientIJ = fullAtmGradientIJ + fullExpGradientIJ;
      const double gradientKI = fullAtmGradientKI + fullExpGradientKI;

      // TODO: guard with newton3 flag
      const double allDrTripledGradientJK = 3. / (allDrTripled * drJKSquare);
      const double cosinesGradientJK = (3./4.) * ((numerator / drJKSquare + numeratorKI * numeratorIJ - numeratorJK * numeratorIJ - numeratorJK * numeratorKI) / allDrSquare);
      const double fullAtmGradientJK = _constC * ((1. + cosines) * allDrTripledGradientJK + cosinesGradientJK / allDrTripled);

      double jkSum = 0.;
      for (size_t n = 0; n < sumFactors.size(); ++n) {
        jkSum += sumFactors.at(n) * (2.*n / (3.*drJK) - _alpha);
      }

      const double fullExpGradientJK = expTerm * (-(1. + cosines) * jkSum / drJK + cosinesGradientJK * sum);
      
      const double gradientJK = fullAtmGradientJK + fullExpGradientJK;

      return std::array<double, 3>{gradientIJ, gradientJK, gradientKI};
    }
  
 private:

 /* Triplet Potential Parameters */
  const double _constC = 1.6152500e-3;
  const double _alpha = 1.3783820e1;
  const std::array<double, 6> _A = {
    -0.3081304e8,
    -0.3519442e10,
    0.4928052e11,
    -0.2182411e12,
    0.3430880e12,
    0
  };
  
 /* Pair Potential Parameters */
  const double _constA = 0.3200711798e8;
  const double _alpha1 = -0.2430565544e2;
  const double _alpha2 = -0.1435536209e2;
  const double _alphaInv1 = -0.4532273868e-1;
  const double _constb = 0.2786344368e2;
  const double _constC6 = 0.8992209265;
  const double _constC8 = 0.7316713603e-1;
  const double _constC10 = 0.7835488511e-2;
  const double _constC12 = 1.1043747590e-3;
  const double _constC14 = 2.0486474980e-4;
  const double _constC16 = 5.0017084700e-5;

  const double _constATilde = 0.8268005465e7;
  const double _alphaTilde = 0.1682493666e1;

  const double _minDistance = 1.2047406;

  const std::array<double, 17> _invFactorials = {1.,
                                                 1.,
                                                 0.5,
                                                 1. / 6.,
                                                 1. / 24.,
                                                 1. / 120.,
                                                 1. / 720.,
                                                 1. / 5040.,
                                                 1. / 40320.,
                                                 1. / 362880.,
                                                 1. / 3628800.,
                                                 1. / 39916800.,
                                                 1. / 479001600.,
                                                 1. / 6227020800.,
                                                 1. / 87178291200.,
                                                 1. / 1307674368000.,
                                                 1. / 20922789888000.};
};

}  // namespace mdLib