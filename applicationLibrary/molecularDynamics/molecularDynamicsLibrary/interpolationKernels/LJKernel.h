/**
 * @file LJKernel.h
 *
 * @date 29.04.2025
 * @author Luis Gall
 */

#pragma once

#include "PairwiseKernel.h"

namespace mdLib {

class LJKernel : public PairwiseKernel<LJKernel> {
 public:
  explicit LJKernel(ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary) : PairwiseKernel<LJKernel>() {
    _PPLibrary = &particlePropertiesLibrary;
  };

  double calculateDerivative(double dr) final {
    return 0.;
  }

  double calculate(double dr) final {
    double dr2 = dr * dr;
    double invdr2 = 1. / dr2;
    double lj6 = _PPLibrary->getMixingSigmaSquared(0, 0) * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = _PPLibrary->getMixing24Epsilon(0, 0) * (lj12 + lj12m6) * invdr2;

    return fac;
  }

 private:
  ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;
};

}  // namespace mdLib