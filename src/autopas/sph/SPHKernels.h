/**
 * @file SPHKernels.h
 * @author seckler
 * @date 22.01.18
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Constants.h"

namespace autopas {
namespace sph {
/**
 * class to define a kernel function for SPH simulations
 */
class SPHKernels {
  /// @todo c++20: make constexpr, once getPI is constexpr.
  static inline const double pi{autopas::utils::Math::getPI()};
  static constexpr double kernelSupportRadius = 2.5;
  // const double C_CFL = 0.3;
 public:
  /**
   * Get the kernelSupportRadius
   * @return the kernel support radius
   */
  static inline double getKernelSupportRadius() { return kernelSupportRadius; }

  /**
   * A kernel function for SPH simulations
   *
   * @param dr distance vector
   * @param h relative kernel support radius
   * @return value of the kernel function
   */
  static inline double W(const std::array<double, 3> dr, const double h) {
    const double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
    return W(dr2, h);
  }

  /**
   * A kernel function for sph simulations
   * @param dr2 squared absolute distance
   * @param h relative kernel support radius
   * @return value of the kernel function
   */
  static inline double W(const double dr2, const double h) {
    const double H = kernelSupportRadius * h;
    if (dr2 < H * H) {
      const double s = sqrt(dr2) / H;  // sqrt(dr * dr) / H;
      const double s1 = 1.0 - s;       // (1.0 - s < 0) ? 0 : 1.0 - s;
      const double s2 = std::max(0., 0.5 - s);
      double r_value = (s1 * s1 * s1) - 4.0 * (s2 * s2 * s2);  //  pow(s1, 3) - 4.0 * pow(s2, 3);
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
  static inline std::array<double, 3> gradW(const std::array<double, 3> dr, const double h) {
    const double H = kernelSupportRadius * h;
    const double drabs = autopas::utils::ArrayMath::L2Norm(dr);
    const double s = drabs / H;  // sqrt(dr * dr) / H;
    const double s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
    const double s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
    double r_value = -3.0 * (s1 * s1) + 12.0 * (s2 * s2);  // - 3.0 * pow(s1, 2) + 12.0 * pow(s2, 2);
    // if # of dimension == 3
    r_value *= 16.0 / pi / (H * H * H);
    const double scale = r_value / (drabs * H + 1.0e-6 * h);
    return autopas::utils::ArrayMath::mulScalar(dr, scale);  // dr * r_value / (sqrt(dr * dr) * H + 1.0e-6 * h);
  }
};
}  // namespace sph
}  // namespace autopas
