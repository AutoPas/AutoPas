/**
 * @file CubeUniform.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Object.h"
#include "autopas/utils/ArrayMath.h"
#include "generators/src/UniformGenerator.h"

/**
 * Class describing a cuboid object filled with uniformly randomly distributed particles.
 */
class CubeUniform : public Object {
 public:
  /**
   * Constructor based on the number of particles.
   * @param velocity
   * @param typeId
   * @param numParticles
   * @param boxLength
   * @param bottomLeftCorner
   */
  CubeUniform(const std::array<double, 3> &velocity, const size_t typeId, const size_t numParticles,
              const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner)
      : CubeUniform(velocity, typeId, boxLength, bottomLeftCorner) {
    _numParticles = numParticles;
    const double volume = boxLength[0] * boxLength[1] * boxLength[2];
    _density = static_cast<double>(numParticles) / volume;
  }

  /**
   * Constructor based on particle density.
   * @param velocity
   * @param typeId
   * @param boxLength
   * @param bottomLeftCorner
   * @param density Target particle density.
   */
  CubeUniform(const std::array<double, 3> &velocity, const size_t typeId, const std::array<double, 3> &boxLength,
              const std::array<double, 3> &bottomLeftCorner, const double density)
      : CubeUniform(velocity, typeId, boxLength, bottomLeftCorner) {
    const double volume = boxLength[0] * boxLength[1] * boxLength[2];
    _numParticles = static_cast<size_t>(std::round(density * volume));
    const double newDensity = static_cast<double>(_numParticles) / volume;
    _density = newDensity;
    if (std::abs(newDensity - density) > 1e-10) {
      std::cout << "CubeUniform: The requested density of " << density << " could not be achieved with the given box "
                << "length. Actual density: " << newDensity << "." << std::endl;
    }
  }

  /**
   * Returns the particle density.
   * @return density of particles.
   */
  [[nodiscard]] double getParticleDensity() const { return _density; }

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
    using namespace autopas::utils::ArrayMath::literals;
    return _bottomLeftCorner + _boxLength;
  }

  /**
   * Converts the object to a human readable string.
   * @return human readable string of the uniform cube.
   */
  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "numberOfParticles"
           << ":  " << _numParticles << "\n";
    output << std::setw(_valueOffset) << std::left << "particle-density"
           << ":  " << _density << "\n";
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << "\n";
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << "\n";
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates the particles based on the configuration of the cube object defined in the yaml file.
   * @param particles The container where the generated particles will be stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    // Wrapper so that std::vector can be used as an AutoPas::ParticleContainer
    auto particlesWrapper = autopasTools::PseudoContainer(particles);

    using namespace autopas::utils::ArrayMath::literals;
    const auto boxMax = _bottomLeftCorner + _boxLength;

    // dummy particle used as a template with id of the first newly generated one
    const ParticleType dummyParticle = getDummyParticle(particles.size());

    autopasTools::generators::UniformGenerator::fillWithParticles(particlesWrapper, dummyParticle, _bottomLeftCorner,
                                                                  boxMax, _numParticles);
  }

 private:
  /**
   * Internal constructor.
   * @param velocity
   * @param typeId
   * @param boxLength
   * @param bottomLeftCorner
   */
  CubeUniform(const std::array<double, 3> &velocity, const size_t typeId, const std::array<double, 3> &boxLength,
              const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId), _boxLength(boxLength), _bottomLeftCorner(bottomLeftCorner) {}

  /**
   * The number of particles in the object.
   */
  size_t _numParticles{};

  /**
   * The length of the box in each direction.
   */
  std::array<double, 3> _boxLength;

  /**
   * The Coordinates of the bottom left front corner.
   */
  std::array<double, 3> _bottomLeftCorner;

  /**
   * Target particle density.
   */
  double _density{};
};
