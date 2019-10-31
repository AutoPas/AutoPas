/**
 * @file CubeGrid.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Objects.h"

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
                 * Getter for ParticlesPerDim
                 * @return particlePerDim
                 */
                [[nodiscard]] const std::array<size_t, 3> &
            getParticlesPerDim() const {
    return particlesPerDim;
  }

  /**
   * Getter for ParticleSpacing
   * @return particleSpacing
   */
  [[nodiscard]] double getParticleSpacing() const { return particleSpacing; }

      /**
       * Getter for total number of Particles for object
       * @return particlesTotal
       */
      [[nodiscard]] size_t getParticlesTotal() const override {
    return std::accumulate(std::begin(particlesPerDim), std::end(particlesPerDim), 1, std::multiplies<double>());
  }

  /**
   * Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   */
  const std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  /**
   * Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   */
  const std::array<double, 3> getBoxMax() const override {
    std::array<double, 3> dppD;
    // copy for type conversion
    std::copy(std::begin(particlesPerDim), std::end(particlesPerDim), std::begin(dppD));
    return autopas::ArrayMath::add(bottomLeftCorner, (autopas::ArrayMath::mulScalar(dppD, particleSpacing)));
  }
  /**
   * Prints the Configuration of the current Object
   */
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
