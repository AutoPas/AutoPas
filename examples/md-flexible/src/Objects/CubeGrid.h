/**
 * @file CubeGrid.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include <functional>
#include <numeric>

#include "Objects.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Class describing a regular 3D particle grid object.
 */
class CubeGrid : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * @param particlesPerDim
   * @param particleSpacing
   * @param bottomLeftCorner
   */
  CubeGrid(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
           const std::array<size_t, 3> &particlesPerDim, double particleSpacing,
           const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        particlesPerDim(particlesPerDim),
        particleSpacing(particleSpacing),
        bottomLeftCorner(bottomLeftCorner) {}

  /**
   * Getter for ParticleSpacing
   * @return particleSpacing
   */
  [[nodiscard]] double getParticleSpacing() const { return particleSpacing; }

  /**
   * Getter for ParticlesPerDim
   * @return particlePerDim
   */
  [[nodiscard]] const std::array<size_t, 3> &getParticlesPerDim() const { return particlesPerDim; }

  [[nodiscard]] size_t getParticlesTotal() const override {
    return std::accumulate(std::begin(particlesPerDim), std::end(particlesPerDim), 1, std::multiplies<double>());
  }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    auto particlesPerDimDouble = autopas::utils::ArrayUtils::static_cast_array<double>(particlesPerDim);
    // subtract one because the first particle is at bottomLeftCorner
    auto particlesPerDimSubOne = autopas::utils::ArrayMath::subScalar(particlesPerDimDouble, 1.);
    auto lastParticleRelative = autopas::utils::ArrayMath::mulScalar(particlesPerDimSubOne, particleSpacing);
    auto lastParticleAbsolute = autopas::utils::ArrayMath::add(bottomLeftCorner, lastParticleRelative);

    return lastParticleRelative;
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "particles-per-dimension"
           << ":  " << autopas::utils::ArrayUtils::to_string(particlesPerDim) << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << particleSpacing << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

 private:
  std::array<size_t, 3> particlesPerDim;
  double particleSpacing;
  std::array<double, 3> bottomLeftCorner;
};
