/**
 * @file CubeUniform.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include <random>

#include "Object.h"
#include "autopas/utils/ArrayMath.h"

#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
#include "AbsoluteMultiSiteMoleculeInitializer.h"
#endif

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
  CubeUniform(const std::array<double, 3> &velocity, unsigned long typeId, size_t numParticles,
              const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner)
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
           << ":  " << _numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates the particles based on the configuration of the cube object defined in the yaml file.
   * @param particles The container where the generated particles will be stored.
   */
#if (MD_FLEXIBLE_MODE != MULTISITE) or not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
  void generate(std::vector<ParticleType> &particles, const std::shared_ptr<const ParticlePropertiesLibraryType>& ppl) const override {
#else
  void generate(std::vector<ParticleType> &particles, const std::shared_ptr<const ParticlePropertiesLibraryType>& ppl,
              MoleculeContainer& moleculeContainer) const override{
#endif

    // Set up random number generation
    std::random_device randomDevice;
    std::mt19937 randomNumberEngine(randomDevice());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (unsigned long i = 0; i < _numParticles; ++i) {

      const std::array<double, 3> position{_bottomLeftCorner[0] + distribution(randomNumberEngine) * _boxLength[0],
                     _bottomLeftCorner[1] + distribution(randomNumberEngine) * _boxLength[1],
                     _bottomLeftCorner[2] + distribution(randomNumberEngine) * _boxLength[2]};
#if (MD_FLEXIBLE_MODE != MULTISITE) or not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
      insertMolecule(position, ppl, particles);
#else
      insertMolecule(position, ppl, particles, moleculeContainer);
#endif
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
