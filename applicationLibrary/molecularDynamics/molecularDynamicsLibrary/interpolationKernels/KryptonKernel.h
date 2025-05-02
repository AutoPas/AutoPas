/**
 * @file KryptonKernel.h
 *
 * @date 29.04.2025
 * @author Luis Gall
 */

#pragma once

#include "PairwiseKernel.h"
#include "autopas/utils/ArrayMath.h"

namespace mdLib {

class KryptonKernel : public PairwiseKernel<KryptonKernel> {
 public:
  explicit KryptonKernel()
      : PairwiseKernel<KryptonKernel>(){

        };

  double calculate(double dr) final {
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

    return (dr >= _minDistance ? firstTerm + secondTerm
                               : _constATilde * distInv2 * std::exp(-_alphaTilde * dr) * (distInv + _alphaTilde));
  }

 private:
  const double _constA = 0.3200711798e8;
  const double _alpha1 = -0.2430565544e1;
  const double _alpha2 = -0.1435536209;
  const double _alphaInv1 = -0.4532273868;
  const double _constb = 0.2786344368e1;
  const double _constC6 = 0.8992209265e6;
  const double _constC8 = 0.7316713603e7;
  const double _constC10 = 0.7835488511e8;
  const double _constC12 = 1.1043747590e9;
  const double _constC14 = 2.0486474980e10;
  const double _constC16 = 5.0017084700e11;

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