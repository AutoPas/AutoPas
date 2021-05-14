/**
 * @file CubeUniform.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include <ctime>

#include "autopas/utils/ArrayMath.h"
#include "Object.h"
#include "src/ParticleAttributes.h"

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
        bottomLeftCorner(bottomLeftCorner) {}

  [[nodiscard]] size_t getParticlesTotal() const override { return numParticles; }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength);
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "numberOfParticles"
           << ":  " << numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  void generate(std::vector<ParticleAttributes> particles) const override {
    ParticleAttributes particle = getDummyParticle(particles.size());
		std::srand(std::time(0));
  	for (unsigned long i = 0; i < numParticles; ++i) {
      particle.id++;
    	particle.positionX = static_cast<double>(std::rand()) / RAND_MAX;
    	particle.positionY = static_cast<double>(std::rand()) / RAND_MAX;
    	particle.positionZ = static_cast<double>(std::rand()) / RAND_MAX;
      particles.push_back(particle);
  	}
  }

 private:
  size_t numParticles;
  std::array<double, 3> boxLength;
  std::array<double, 3> bottomLeftCorner;
};
