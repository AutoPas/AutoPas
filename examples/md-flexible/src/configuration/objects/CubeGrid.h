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
#include "generators/src/GridGenerator.h"

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
   * @param density Target particle density (optional, computed from spacing if 0).
   * @param centered If true, the cell offset moved inward by 1/2 * spacing.
   */
  CubeGrid(const std::array<double, 3> &velocity, unsigned long typeId, const std::array<size_t, 3> &particlesPerDim,
           double particleSpacing, const std::array<double, 3> &bottomLeftCorner, const double density = 0.0,
           const bool centered = true)
      : Object(velocity, typeId),
        _particlesPerDim(particlesPerDim),
        _particleSpacing(particleSpacing),
        _bottomLeftCorner(bottomLeftCorner),
        _density(density),
      _centered(centered) {
    if (_particleSpacing <= 0.0 and _density > 0.0) {
      _particleSpacing = std::cbrt(1.0 / _density);
    } else if (_particleSpacing > 0.0 and _density <= 0.0) {
      _density = 1.0 / (_particleSpacing * _particleSpacing * _particleSpacing);
    }
  }

  /**
   * Returns the particle spacing.
   * @return spacing between particles.
   */
  [[nodiscard]] double getParticleSpacing() const override { return _particleSpacing; }

  /**
   * Returns the particle density.
   * @return density of particles.
   */
  [[nodiscard]] double getParticleDensity() const { return _density; }

  /**
   * Getter for ParticlesPerDim
   * @return particlePerDim
   */
  [[nodiscard]] const std::array<size_t, 3> &getParticlesPerDim() const { return _particlesPerDim; }

  /**
   * Returns the total number of particles which will be / have been generated.
   * @return number of generated particles.
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    return std::accumulate(std::begin(_particlesPerDim), std::end(_particlesPerDim), 1, std::multiplies<double>());
  }

  /**
   * Returns the coordinates of the bottom left front corner.
   * @return bottom left front corner.
   */
  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  /**
   * Returns the coordinates of the top right back corner.
   * @return top right back corner.
   */
  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    using namespace autopas::utils::ArrayMath::literals;

    const auto particlesPerDimDouble = autopas::utils::ArrayUtils::static_cast_copy_array<double>(_particlesPerDim);
    const auto totalLengthRelative = particlesPerDimDouble * _particleSpacing;
    return _bottomLeftCorner + totalLengthRelative;
  }

  /**
   * Turns the cube grid object into a human-readable string.
   * @returns human-readable string of a cube grid object.
   */
  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "particles-per-dimension"
           << ":  " << autopas::utils::ArrayUtils::to_string(_particlesPerDim) << "\n";
    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << _particleSpacing << "\n";
    output << std::setw(_valueOffset) << std::left << "particle-density"
           << ":  " << _density << "\n";
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << "\n";
    output << std::setw(_valueOffset) << std::left << "centered"
           << ":  " << std::to_string(_centered) << "\n";
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates the particles based on the configuration of the CubeGrid object provided in the yaml file.
   * @param particles The container in which the generated particles get stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    using namespace autopas::utils::ArrayMath::literals;

    // Wrapper so that std::vector can be used as an AutoPas::ParticleContainer
    auto particlesWrapper = autopasTools::PseudoContainer(particles);

    // dummy particle used as a template with id of the first newly generated one
    const ParticleType dummyParticle = getDummyParticle(particles.size());

    std::array<double, 3> offset = _bottomLeftCorner;
    if (_centered) {
      offset += 0.5 * _particleSpacing;
    }

    autopasTools::generators::GridGenerator::fillWithParticles(particlesWrapper, _particlesPerDim, dummyParticle,
                                                               {_particleSpacing, _particleSpacing, _particleSpacing},
                                                               offset);
  }

 private:
  /**
   * Defines how many particles will be created in each dimension.
   */
  std::array<size_t, 3> _particlesPerDim;

  /**
   * Defines the amount of space between particles.
   */
  double _particleSpacing;

  /**
   * Stores the coordinates of the bottom left front corner.
   */
  std::array<double, 3> _bottomLeftCorner;

  /**
   * Target particle density.
   */
  double _density{0.0};

  /**
   *
   */
bool _centered{true};
};
