/**
 * @file SphereGenerator.h
 * @author N. Fottner
 * @date 29/10/19
 */

#pragma once

#include <array>
#include <functional>

#include "autopas/utils/ArrayMath.h"

/**
 * Generator for a regular 3D spherical particle grid. Works by starting with a uniform grid and then removing
 * everything which doesn't fit into the sphere.
 */
namespace autopas::generators::SphereGenerator {

/**
 * Call f for every point on the sphere where a particle should be.
 * @param center coordinates of the sphere's center.
 * @param radius radius of the sphere in number of particles.
 * @param particleSpacing The amount of space between each particle.
 * @param f Function called for every point.
 */
inline void iteratePositions(const std::array<double, 3> &center, int radius, double particleSpacing,
                             const std::function<void(std::array<double, 3>)> &f) {
  using namespace autopas::utils::ArrayMath::literals;

  // generate regular grid for 1/8th of the sphere
  for (int z = 0; z <= radius; ++z) {
    for (int y = 0; y <= radius; ++y) {
      for (int x = 0; x <= radius; ++x) {
        // position relative to the center
        const std::array<double, 3> relativePos = {(double)x, (double)y, (double)z};
        // mirror to rest of sphere
        for (int i = -1; i <= 1; i += 2) {
          for (int k = -1; k <= 1; k += 2) {
            for (int l = -1; l <= 1; l += 2) {
              const std::array<double, 3> mirrorMultipliers = {(double)i, (double)k, (double)l};
              // position mirrored, scaled and absolute
              const std::array<double, 3> posVector = center + ((relativePos * mirrorMultipliers) * particleSpacing);

              double distFromCentersSquare = autopas::utils::ArrayMath::dot(posVector - center, posVector - center);
              const auto r = (radius + 1) * particleSpacing;
              const auto rSquare = r * r;
              // since the loops create a cubic grid only apply f for positions inside the sphere
              if (distFromCentersSquare <= rSquare) {
                f(posVector);
              }
              // avoid duplicates
              if (z == 0) break;
            }
            if (y == 0) break;
          }
          if (x == 0) break;
        }
      }
    }
  }
}

}  // namespace autopas::generators::SphereGenerator