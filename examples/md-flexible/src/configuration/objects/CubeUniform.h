/**
 * @file CubeUniform.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Object.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/UniformGenerator.h"

/**
 * Class describing an cuboid object filled with uniformly randomly distributed particles.
 */
class CubeUniform : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param numParticles
   * @param boxLength
   * @param bottomLeftCorner
   */
  CubeUniform(const std::array<CalcPrecision, 3> &velocity, unsigned long typeId, size_t numParticles,
              const std::array<CalcPrecision, 3> &boxLength, const std::array<CalcPrecision, 3> &bottomLeftCorner)
      : Object(velocity, typeId),
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
  [[nodiscard]] std::array<CalcPrecision, 3> getBoxMin() const override { return _bottomLeftCorner; }

  /**
   * Returns the coordinates of the top right back corner.
   * @return top right back corner of the cube.
   */
  [[nodiscard]] std::array<CalcPrecision, 3> getBoxMax() const override {
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
   * The number of particles in the object.
   */
  size_t _numParticles;

  /**
   * The lenght of the box in each direction.
   */
  std::array<CalcPrecision, 3> _boxLength;

  /**
   * The Coordinates of the bottom left front corner.
   */
  std::array<CalcPrecision, 3> _bottomLeftCorner;
};
