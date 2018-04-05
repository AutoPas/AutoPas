//
// Created by seckler on 22.01.18.
//

#ifndef AUTOPAS_SPHKERNELS_H
#define AUTOPAS_SPHKERNELS_H

#include "utils/arrayMath.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace autopas {
namespace sph {
/**
 * class to define a kernel function for SPH simulations
 */
class SPHKernels {
  static constexpr double pi = M_PI;  // atan(1.0) * 4.0;
  static constexpr double kernelSupportRadius = 2.5;
  // const double C_CFL = 0.3;
 public:

  static inline double getKernelSupportRadius(){
    return kernelSupportRadius;
  }

  /**
   * A kernel function for SPH simulations
   *
   * @param dr distance vector
   * @param h relative kernel support radius
   * @return value of the kernel function
   */
  static inline double W(const std::array<double, 3> dr, const double h) {
    const double H = kernelSupportRadius * h;
    const double dr2 = autopas::arrayMath::dot(dr, dr);
    if (dr2 < H * H) {
      const double s = sqrt(dr2) / H;  // sqrt(dr * dr) / H;
      const double s1 = 1.0 - s;       // (1.0 - s < 0) ? 0 : 1.0 - s;
      const double s2 = std::max(0., 0.5 - s);
      double r_value = (s1 * s1 * s1) -
                       4.0 * (s2 * s2 * s2);  //  pow(s1, 3) - 4.0 * pow(s2, 3);
      // if # of dimension == 3
      r_value *= 16.0 / pi / (H * H * H);
      return r_value;
    } else {
      return 0.;
    }
  }

  /**
   * returns the flops for one full calculation of the kernel
   * @return flops for one full calculation of the kernel
   */
  static unsigned long getFlopsW();

  /**
   * gradient of the kernel function W
   * @param dr distance vector
   * @param h kernel support radius
   * @return gradient of W evaluated at dr and h
   */
  static inline std::array<double, 3> gradW(const std::array<double, 3> dr,
                                            const double h) {
    const double H = kernelSupportRadius * h;
    const double drabs = sqrt(autopas::arrayMath::dot(dr, dr));
    const double s = drabs / H;  // sqrt(dr * dr) / H;
    const double s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
    const double s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
    double r_value =
        -3.0 * (s1 * s1) +
        12.0 * (s2 * s2);  // - 3.0 * pow(s1, 2) + 12.0 * pow(s2, 2);
    // if # of dimension == 3
    r_value *= 16.0 / pi / (H * H * H);
    const double scale = r_value / (drabs * H + 1.0e-6 * h);
    return autopas::arrayMath::mulScalar(
        dr, scale);  // dr * r_value / (sqrt(dr * dr) * H + 1.0e-6 * h);
  }
};
}  // namespace sph
}  // namespace autopas
#endif  // AUTOPAS_SPHKERNELS_H
