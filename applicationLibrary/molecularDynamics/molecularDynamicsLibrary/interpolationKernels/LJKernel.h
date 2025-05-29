/**
 * @file LJKernel.h
 *
 * @date 29.04.2025
 * @author Luis Gall
 */

#pragma once

#include "Kernel.h"

namespace mdLib {

class LJKernel : public Kernel<LJKernel> {
 public:
  explicit LJKernel(ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary) : Kernel<LJKernel>() {
    _PPLibrary = &particlePropertiesLibrary;
  };

  double calculatePairDerivative(double dr) final {
    return 0.;
  }

  double calculateTriplet(double dr1, double dr2, double dr3) final {
    return 0.;
  }

  double calculatePair(double dr) final {
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