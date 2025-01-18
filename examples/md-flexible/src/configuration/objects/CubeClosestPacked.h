/**
 * @file CubeClosestPacked.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <cmath>

#include "Object.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/PseudoContainer.h"
#include "autopasTools/generators/ClosestPackingGenerator.h"

/**
 * Class describing a cube of hexagonally closest packed particles.
 */
class CubeClosestPacked : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param particleSpacing distance between all neighboring particles
   * @param boxLength
   * @param bottomLeftCorner
   */
  CubeClosestPacked(const std::array<CalcType, 3> &velocity, unsigned long typeId, CalcType particleSpacing,
                    const std::array<CalcType, 3> &boxLength, const std::array<CalcType, 3> &bottomLeftCorner)
      : Object(velocity, typeId),
        _boxLength(boxLength),
        _particleSpacing(particleSpacing),
        _bottomLeftCorner(bottomLeftCorner),
        _topRightCorner(autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength)),
        _xOffset(particleSpacing * 1. / 2.),
        _yOffset(particleSpacing * sqrt(1. / 12.)) {}

  [[nodiscard]] CalcType getParticleSpacing() const override { return _particleSpacing; }

  /**
   * Returns the total amount of particles which will be / have been generated.
   * @return number of generated particles.
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    // Number of particles in the first row.
    const size_t xNumRow = std::ceil(_boxLength[0] / _particleSpacing);
    // True if the total number of x-positions is odd.
    const bool xOdd = static_cast<int>(std::ceil(_boxLength[0] / _xOffset)) % 2 == 1;

    // Distance between layers and rows
    const auto spacingLayer = _particleSpacing * sqrt(2. / 3.);
    const auto spacingRow = _particleSpacing * sqrt(3. / 4.);

    // Number of rows in an even layer.
    const size_t yNumEven = std::ceil(_boxLength[1] / spacingRow);
    // Number of rows in an odd layer.
    const size_t yNumOdd = std::ceil((_boxLength[1] - _yOffset) / spacingRow);

    // Number of particles in an even layer.
    const size_t evenLayer = xNumRow * yNumEven - std::floor(xOdd * yNumEven * 0.5);
    // Number of particles in an odd layer.
    const size_t oddLayer = xNumRow * yNumOdd - std::ceil(xOdd * yNumOdd * 0.5);

    // Total number of layers.
    const CalcType numLayers = std::ceil(_boxLength[2] / spacingLayer);
    // Add up all even and odd layers.
    return evenLayer * std::ceil(numLayers / 2.) + oddLayer * std::floor(numLayers / 2.);
  }

  [[nodiscard]] std::array<CalcType, 3> getBoxMin() const override { return _bottomLeftCorner; }

  [[nodiscard]] std::array<CalcType, 3> getBoxMax() const override {
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
           << ":  " << _particleSpacing << "\n";
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << "\n";
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << "\n";
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates particles based on the parameters provided to the CubeClosestPacked Object in the configuration file.
   * @param particles: The container, where the new particles get stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    // Wrapper so that std::vector can be used as an AutoPas::ParticleContainer
    auto particlesWrapper = autopasTools::PseudoContainer(particles);

    // dummy particle used as a template with id of the first newly generated one
    const ParticleType dummyParticle = getDummyParticle(particles.size());

    autopasTools::generators::ClosestPackingGenerator::fillWithParticles(
        particlesWrapper, _bottomLeftCorner, _topRightCorner, dummyParticle, _particleSpacing);
  }

 private:
  /**
   * The distance between the particles.
   */
  CalcType _particleSpacing;

  /**
   * Extend of the box in each dimension.
   */
  std::array<CalcType, 3> _boxLength;

  /**
   * Minimum box coordinates.
   */
  std::array<CalcType, 3> _bottomLeftCorner;

  /**
   * Maximum box coordinates
   */
  std::array<CalcType, 3> _topRightCorner;

  /**
   * Shorter part of the bisectrix when split at the intersection of all bisectrices.
   */
  CalcType _xOffset;

  /**
   * Shorter part of the bisectrix when split at the intersection of all bisectrices.
   */
  CalcType _yOffset;
};
