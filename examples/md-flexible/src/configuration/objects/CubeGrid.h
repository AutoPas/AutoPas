/**
 * @file CubeGrid.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include <functional>
#include <numeric>

#include "src/ParticleAttributes.h"
#include "autopas/utils/ArrayMath.h"
#include "Object.h"

/**
 * Class describing a regular 3D particle grid object.
 */
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
   * @param bottomLeftCorner
   */
  CubeGrid(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
           const std::array<size_t, 3> &particlesPerDim, double particleSpacing,
           const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        _particlesPerDim(particlesPerDim),
        _particleSpacing(particleSpacing),
        _bottomLeftCorner(bottomLeftCorner) {}

  [[nodiscard]] double getParticleSpacing() const override { return _particleSpacing; }

  /**
   * Getter for ParticlesPerDim
   * @return particlePerDim
   */
  [[nodiscard]] const std::array<size_t, 3> &getParticlesPerDim() const { return _particlesPerDim; }

  [[nodiscard]] size_t getParticlesTotal() const override {
    return std::accumulate(std::begin(_particlesPerDim), std::end(_particlesPerDim), 1, std::multiplies<double>());
  }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    auto particlesPerDimDouble = autopas::utils::ArrayUtils::static_cast_array<double>(_particlesPerDim);
    // subtract one because the first particle is at bottomLeftCorner
    auto particlesPerDimSubOne = autopas::utils::ArrayMath::subScalar(particlesPerDimDouble, 1.);
    auto lastParticleRelative = autopas::utils::ArrayMath::mulScalar(particlesPerDimSubOne, _particleSpacing);
    auto lastParticleAbsolute = autopas::utils::ArrayMath::add(_bottomLeftCorner, lastParticleRelative);

    return lastParticleAbsolute;
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "particles-per-dimension"
           << ":  " << autopas::utils::ArrayUtils::to_string(_particlesPerDim) << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << _particleSpacing << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  void generate(std::vector<ParticleAttributes> &particles) const override {
    ParticleAttributes particle = getDummyParticle(particles.size());
		std::array<double, 3> offset = {0.0, 0.0, 0.0};

		for (unsigned long z = 0; z < _particlesPerDim[2]; ++z) {
			for (unsigned long y = 0; y < _particlesPerDim[1]; ++y) {
				for (unsigned long x = 0; x < _particlesPerDim[0]; ++x) {
          particle.id++;
    			particle.position[0] = (static_cast<double>(x) * _particleSpacing + offset[0]);
    			particle.position[1] = (static_cast<double>(y) * _particleSpacing + offset[1]);
    			particle.position[2] = (static_cast<double>(z) * _particleSpacing + offset[2]);

          particles.push_back(particle);
				}
			}
		}
  }

 private:
  std::array<size_t, 3> _particlesPerDim;
  double _particleSpacing;
  std::array<double, 3> _bottomLeftCorner;
};
