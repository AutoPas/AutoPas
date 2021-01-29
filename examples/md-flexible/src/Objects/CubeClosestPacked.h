/**
 * @file CubeClosestPacked.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <math.h>

#include <functional>

#include "Objects.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Class describing a cube of hexagonally closest packed particles.
 */
class CubeClosestPacked : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * @param particleSpacing distance between all neighboring particles
   * @param boxLength
   * @param bottomLeftCorner
   */
  CubeClosestPacked(const std::array<double, 3> &velocity, uint64_t typeId, double epsilon, double sigma,
                    double mass, double particleSpacing, const std::array<double, 3> &boxLength,
                    const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        boxLength(boxLength),
        particleSpacing(particleSpacing),
        bottomLeftCorner(bottomLeftCorner) {}

  [[nodiscard]] double getParticleSpacing() const override { return particleSpacing; }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength);
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << particleSpacing << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

 private:
  double particleSpacing;
  std::array<double, 3> boxLength;
  std::array<double, 3> bottomLeftCorner;
};
