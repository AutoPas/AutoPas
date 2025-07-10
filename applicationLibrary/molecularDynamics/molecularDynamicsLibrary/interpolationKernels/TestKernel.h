/**
 * @file TestKernel.h
 *
 * @date 25.06.2025
 * @author Luis Gall
 */

#pragma once

#include "Kernel.h"
#include "autopas/utils/ArrayMath.h"

namespace mdLib {

class TestKernel : public Kernel<TestKernel> {
 public:
  explicit TestKernel()
      : Kernel<TestKernel>(){

        };

  double calculatePairDerivative(double dr) final {
    return 0.;
  }

  double calculatePair(double dr) final {
    return 0.;
  }

  double calculateTriplet(double dr1, double dr2, double dr3) {
    return dr1 * dr1 + dr2 * dr2 + dr3 * dr3;
  }

  std::array<double, 3> calculateTripletDerivative(double dr1, double dr2, double dr3) {
    return std::array{2.*dr1 + 2.*dr2 + 2.*dr3, 0., 0.};
  }

};

}  // namespace mdLib