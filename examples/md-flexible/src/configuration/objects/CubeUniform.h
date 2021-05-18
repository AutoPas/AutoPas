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
   * @param _numParticles
   * @param _boxLength
   * @param _bottomLeftCorner
   */
  CubeUniform(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
              size_t _numParticles, const std::array<double, 3> &_boxLength,
              const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        _numParticles(_numParticles),
        _boxLength(_boxLength),
        _bottomLeftCorner(bottomLeftCorner) {}

  [[nodiscard]] size_t getParticlesTotal() const override { return _numParticles; }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(_bottomLeftCorner, _boxLength);
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "numberOfParticles"
           << ":  " << _numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "_bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  void generate(std::vector<ParticleAttributes> &particles) const override {
    ParticleAttributes particle = getDummyParticle(particles.size());
		std::srand(std::time(0));
  	for (unsigned long i = 0; i < _numParticles; ++i) {
      particle.id++;
    	particle.positionX = _bottomLeftCorner[0] + (static_cast<double>(std::rand()) / RAND_MAX) * _boxLength[0];
    	particle.positionY = _bottomLeftCorner[1] + (static_cast<double>(std::rand()) / RAND_MAX) * _boxLength[1];
    	particle.positionZ = _bottomLeftCorner[2] + (static_cast<double>(std::rand()) / RAND_MAX) * _boxLength[2];
      particles.push_back(particle);
  	}
  }

 private:
  size_t _numParticles;
  std::array<double, 3> _boxLength;
  std::array<double, 3> _bottomLeftCorner;
};
