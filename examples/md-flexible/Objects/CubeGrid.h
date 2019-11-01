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
   * @param particlesTotal
   * @param bottomLeftCorner
   */
  CubeGrid(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
           const std::array<size_t, 3> &particlesPerDim, double particleSpacing,
           const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        particlesPerDim(particlesPerDim),
        particleSpacing(particleSpacing),
        bottomLeftCorner(bottomLeftCorner){}

            /**
             * Getter for ParticleSpacing
             * @return particleSpacing
             */
            [[nodiscard]] double getParticleSpacing() const {
    return particleSpacing;
  }

  /**
   * Getter for ParticlesPerDim
   * @return particlePerDim
   */
  [[nodiscard]] const std::array<size_t, 3> &getParticlesPerDim() const { return particlesPerDim; }

      [[nodiscard]] size_t getParticlesTotal() const override {
    return std::accumulate(std::begin(particlesPerDim), std::end(particlesPerDim), 1, std::multiplies<double>());
  }

  const std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  const std::array<double, 3> getBoxMax() const override {
    std::array<double, 3> dppD;
    // copy for type conversion
    std::copy(std::begin(particlesPerDim), std::end(particlesPerDim), std::begin(dppD));
    return autopas::ArrayMath::add(bottomLeftCorner, (autopas::ArrayMath::mulScalar(dppD, particleSpacing)));
  }

  std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "Particles per dimension"
           << ":  " << autopas::ArrayUtils::to_string(particlesPerDim) << std::endl;
    output << std::setw(_valueOffset) << std::left << "Particle spacing"
           << ":  " << particleSpacing << std::endl;
    output << std::setw(_valueOffset) << std::left << "Number of Particles"
           << ":  " << (particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]) << std::endl;
    output << Object::to_string();
    return output.str();
  }

 private:
  std::array<size_t, 3> particlesPerDim;
  double particleSpacing;
  std::array<double, 3> bottomLeftCorner;
};
