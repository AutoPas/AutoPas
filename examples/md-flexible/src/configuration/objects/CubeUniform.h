/**
 * @file CubeUniform.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include <ctime>

#include "Object.h"
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
        _numParticles(numParticles),
        _boxLength(boxLength),
        _bottomLeftCorner(bottomLeftCorner) {}

  /**
   * Returns the total amount of particles which will be / have been generated.
   * @return total amount of particles.
   */
  [[nodiscard]] size_t getParticlesTotal() const override { return _numParticles; }

  /**
   * Returns the coordinates of the bottom left front corner.
   * @return bottom left front corner of the cube.
   */
  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  /**
   * Returns the coordinates of the top right back corner.
   * @return top right back corner of the cube.
   */
  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(_bottomLeftCorner, _boxLength);
  }

  /**
   * Converts the object to a human readable string.
   * @return human readable string of the uniform cube.
   */
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

  /**
   * Generates the particles based on the configuration of the cube object defined in the yaml file.
   * @param particles The container where the generated particles will be stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    ParticleType particle = getDummyParticle(particles.size());
    std::srand(std::time(0));
    for (unsigned long i = 0; i < _numParticles; ++i) {
      particle.setR({_bottomLeftCorner[0] + (static_cast<double>(std::rand()) / RAND_MAX) * _boxLength[0],
                     _bottomLeftCorner[1] + (static_cast<double>(std::rand()) / RAND_MAX) * _boxLength[1],
                     _bottomLeftCorner[2] + (static_cast<double>(std::rand()) / RAND_MAX) * _boxLength[2]});
      particles.push_back(particle);
      particle.setID(particle.getID() + 1);
    }
  }

 private:
  /**
   * The number of particles in the object.
   */
  size_t _numParticles;

  /**
   * The lenght of the box in each direction.
   */
  std::array<double, 3> _boxLength;

  /**
   * The Coordinates of the bottom left front corner.
   */
  std::array<double, 3> _bottomLeftCorner;
};
