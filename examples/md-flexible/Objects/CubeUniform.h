/**
 * @file CubeUniform.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Objects.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Class describing an cuboid object filled with uniformly randomly distributed particles.
 */
class CubeUniform : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * @param numParticles
   * @param boxLength
   * @param bottomLeftCorner
   */
  CubeUniform(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
              size_t numParticles, const std::array<double, 3> &boxLength,
              const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        numParticles(numParticles),
        boxLength(boxLength),
        bottomLeftCorner(bottomLeftCorner){}

            [[nodiscard]] size_t getParticlesTotal() const override {
    return numParticles;
  }

  const std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  const std::array<double, 3> getBoxMax() const override {
    return autopas::ArrayMath::add(bottomLeftCorner, boxLength);
  }

  std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "Center"
           << ":  " << autopas::ArrayUtils::to_string(bottomLeftCorner) << std::endl;
    output << std::setw(_valueOffset) << std::left << "NumberOfParticles"
           << ":  " << numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "BoxLength"
           << ":  " << autopas::ArrayUtils::to_string(boxLength) << std::endl;
    output << Object::to_string();
    return output.str();
  }

 private:
  size_t numParticles;
  std::array<double, 3> boxLength;
  std::array<double, 3> bottomLeftCorner;
};