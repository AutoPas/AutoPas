/**
 * @file CubeClosestPacked.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <cmath>
#include <functional>

#include "Object.h"
#include "autopas/utils/ArrayMath.h"

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
  CubeClosestPacked(const std::array<double, 3> &velocity, unsigned long typeId, double particleSpacing,
                    const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId),
        _boxLength(boxLength),
        _particleSpacing(particleSpacing),
        _bottomLeftCorner(bottomLeftCorner),
        _topRightCorner(autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength)),
        _spacingLayer(particleSpacing * sqrt(2. / 3.)),
        _spacingRow(particleSpacing * sqrt(3. / 4.)),
        _xOffset(particleSpacing * 1. / 2.),
        _yOffset(particleSpacing * sqrt(1. / 12.)) {}

  [[nodiscard]] double getParticleSpacing() const override { return _particleSpacing; }

  /**
   * Returns the total amount of particles which will be / have been generated.
   * @return number of generated particles.
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    // Number of particles in the first row.
    const size_t xNumRow = std::ceil(_boxLength[0] / _particleSpacing);
    // True if the total number of x-positions is odd.
    const bool xOdd = static_cast<int>(std::ceil(_boxLength[0] / _xOffset)) % 2 == 1;

    // Number of rows in an even layer.
    const size_t yNumEven = std::ceil(_boxLength[1] / _spacingRow);
    // Number of rows in an odd layer.
    const size_t yNumOdd = std::ceil((_boxLength[1] - _yOffset) / _spacingRow);

    // Number of particles in an even layer.
    const size_t evenLayer = xNumRow * yNumEven - std::floor(xOdd * yNumEven * 0.5);
    // Number of particles in an odd layer.
    const size_t oddLayer = xNumRow * yNumOdd - std::ceil(xOdd * yNumOdd * 0.5);

    // Total number of layers.
    const double numLayers = std::ceil(_boxLength[2] / _spacingLayer);
    // Add up all even and odd layers.
    return evenLayer * std::ceil(numLayers / 2.) + oddLayer * std::floor(numLayers / 2.);
  }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    using namespace autopas::utils::ArrayMath::literals;
    return _bottomLeftCorner + _boxLength;
  }

  /**
   * Converts the object to a human readable string
   * @return the generated string
   */
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

  /**
   * Generates particles based on the parameters provided to the CubeClosestPacked Object in the configuration file.
   * @param particles: The container, where the new particles get stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    ParticleType particle = getDummyParticle(particles.size());

    bool evenLayer = true;

    for (double z = _bottomLeftCorner[2]; z < _topRightCorner[2]; z += _spacingLayer) {
      double starty = evenLayer ? _bottomLeftCorner[1] : _bottomLeftCorner[1] + _yOffset;
      bool evenRow = evenLayer;  // To ensure layers are alternating as for hexagonal close packed.
      for (double y = starty; y < _topRightCorner[1]; y += _spacingRow) {
        double startx = evenRow ? _bottomLeftCorner[0] : _bottomLeftCorner[0] + _xOffset;
        for (double x = startx; x < _topRightCorner[0]; x += _particleSpacing) {
          particle.setR({x, y, z});
          particles.push_back(particle);

          particle.setID(particle.getID() + 1);
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
   * Extend of the box in each dimension.
   */
  std::array<double, 3> _boxLength;

  /**
   * Minimum box coordinates.
   */
  std::array<double, 3> _bottomLeftCorner;

  /**
   * Maximum box coordinates
   */
  std::array<double, 3> _topRightCorner;

  /**
   * Spacing in y direction when only moving 60° on the unit circle. Or the height in an equilateral triangle.
   */
  double _spacingRow;

  /**
   * Spacing in z direction. Height in an equilateral tetraeder.
   */
  double _spacingLayer;

  /**
   * Shorter part of the bisectrix when split at the intersection of all bisectrices.
   */
  double _xOffset;

  /**
   * Shorter part of the bisectrix when split at the intersection of all bisectrices.
   */
  double _yOffset;
};
