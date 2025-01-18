/**
 * @file CubeGrid.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include <functional>
#include <numeric>

#include "Object.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/GridGenerator.h"

/**
 * Class describing a regular 3D particle grid object.
 */
class CubeGrid : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param particlesPerDim
   * @param particleSpacing
   * @param bottomLeftCorner
   */
  CubeGrid(const std::array<CalcType, 3> &velocity, unsigned long typeId, const std::array<size_t, 3> &particlesPerDim,
           CalcType particleSpacing, const std::array<CalcType, 3> &bottomLeftCorner)
      : Object(velocity, typeId),
        _particlesPerDim(particlesPerDim),
        _particleSpacing(particleSpacing),
        _bottomLeftCorner(bottomLeftCorner) {}

  /**
   * Returns the particle spacing.
   * @return spacing between particles.
   */
  [[nodiscard]] CalcType getParticleSpacing() const override { return _particleSpacing; }

  /**
   * Getter for ParticlesPerDim
   * @return particlePerDim
   */
  [[nodiscard]] const std::array<size_t, 3> &getParticlesPerDim() const { return _particlesPerDim; }

  /**
   * Returns the total amount of particles which will be / have been generated.
   * @return number of generated particles.
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    return std::accumulate(std::begin(_particlesPerDim), std::end(_particlesPerDim), 1, std::multiplies<CalcType>());
  }

  /**
   * Returns the coordinates of the bottom left front corner.
   * @return bottom left front corner.
   */
  [[nodiscard]] std::array<CalcType, 3> getBoxMin() const override { return _bottomLeftCorner; }

  /**
   * Returns the coordinates of the top right back corner.
   * @return top right back corner.
   */
  [[nodiscard]] std::array<CalcType, 3> getBoxMax() const override {
    using namespace autopas::utils::ArrayMath::literals;

    const auto particlesPerDimCalcType = autopas::utils::ArrayUtils::static_cast_copy_array<CalcType>(_particlesPerDim);
    // subtract one because the first particle is at bottomLeftCorner
    const auto particlesPerDimSubOne = particlesPerDimCalcType - static_cast<CalcType>(1.);
    const auto lastParticleRelative = particlesPerDimSubOne * _particleSpacing;
    auto lastParticleAbsolute = _bottomLeftCorner + lastParticleRelative;

    return lastParticleAbsolute;
  }

  /**
   * Turns the cube grid object into a human readable string.
   * @returns human readable string of cube grid object.
   */
  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "particles-per-dimension"
           << ":  " << autopas::utils::ArrayUtils::to_string(_particlesPerDim) << "\n";
    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << _particleSpacing << "\n";
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << "\n";
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates the particles based on the configuration of the CubeGrid object provided in the yaml file.
   * @param particles The container in which the generated particles get stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    // Wrapper so that std::vector can be used as an AutoPas::ParticleContainer
    auto particlesWrapper = autopasTools::PseudoContainer(particles);

    // dummy particle used as a template with id of the first newly generated one
    const ParticleType dummyParticle = getDummyParticle(particles.size());

    autopasTools::generators::GridGenerator::fillWithParticles(particlesWrapper, _particlesPerDim, dummyParticle,
                                                               {_particleSpacing, _particleSpacing, _particleSpacing},
                                                               _bottomLeftCorner);
  }

 private:
  /**
   * Defines how many particles will be created in each dimension.
   */
  std::array<size_t, 3> _particlesPerDim;

  /**
   * Defines the amount of space between particles.
   */
  CalcType _particleSpacing;

  /**
   * Stores the coordinates of the bottom left front corner.
   */
  std::array<CalcType, 3> _bottomLeftCorner;
};
