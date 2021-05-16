/**
 * @file CubeClosestPacked.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <math.h>
#include <functional>

#include "src/ParticleAttributes.h"
#include "autopas/utils/ArrayMath.h"
#include "Object.h"

/**
 * Class describing a cube of hexagonally closest packed particles.
 */
class CubeClosestPacked : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * @param particleSpacing distance between all neighboring particles
   * @param boxLength
   * @param bottomLeftCorner
   */
  CubeClosestPacked(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma,
                    double mass, double particleSpacing, const std::array<double, 3> &boxLength,
                    const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        _boxLength(boxLength),
        _particleSpacing(particleSpacing),
        _bottomLeftCorner(bottomLeftCorner) {}

  [[nodiscard]] double getParticleSpacing() const override { return _particleSpacing; }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(_bottomLeftCorner, _boxLength);
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << _particleSpacing << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  void generate(std::vector<ParticleAttributes> &particles) const override {
    ParticleAttributes particle = getDummyParticle(particles.size());

    const double spacingRow = _particleSpacing * sqrt(3. / 4.);
    const double spacingLayer = _particleSpacing * sqrt(2. / 3.);
    const double xOffset = _particleSpacing * 1. / 2.;
    const double yOffset = _particleSpacing * sqrt(1. / 12.);

    bool evenLayer = true;
    bool evenRow = true;

    for (double z = 0; z < _boxLength[2]; z += spacingLayer) {

      double starty = evenLayer ? 0 : yOffset;
      for (double y = starty; y < _boxLength[1]; y += spacingRow) {

        double startx = evenRow ? 0 : xOffset;
        for (double x = startx; x < _boxLength[0]; x += _particleSpacing) {
    			particle.positionX = x;
    			particle.positionY = y;
    			particle.positionZ = z;
          particles.push_back(particle);

          particle.id++;
        }
        evenRow = not evenRow;
      }
      evenLayer = not evenLayer;
    }
	}

 private:
	/**
	* The distance between the particles.
	*/
  double _particleSpacing;

	/**
	* Extend of the box in each dimension>
	*/
  std::array<double, 3> _boxLength;

	/**
	* Location of the boxe's bottom left corner.
	*/
  std::array<double, 3> _bottomLeftCorner;
};
