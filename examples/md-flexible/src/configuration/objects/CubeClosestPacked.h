/**
 * @file CubeClosestPacked.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <cmath>

#include "Object.h"
#include "autopas/utils/ArrayMath.h"
#include "generators/src/FCCGenerator.h"
#include "generators/src/HCPGenerator.h"
#include "generators/src/PseudoContainer.h"

/**
 * Class describing a cube of closest packed particles (FCC or HCP).
 */
class CubeClosestPacked : public Object {
 public:
  /**
   * Structure type of closest packing.
   */
  enum Structure { fcc, hcp };

  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param particleSpacing distance between all neighboring particles
   * @param boxLength
   * @param bottomLeftCorner
   * @param density Target particle density (optional, computed from spacing if 0).
   * @param structure fcc (default) or hcp.
   * @param centered If true, the cell offset moved inward by 1/4 * lattice constant.
   */
  CubeClosestPacked(const std::array<double, 3> &velocity, unsigned long typeId, double particleSpacing,
                    const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner,
                    double density = 0.0, Structure structure = Structure::fcc, const bool centered = true)
      : Object(velocity, typeId),
        _boxLength(boxLength),
        _particleSpacing(particleSpacing),
        _bottomLeftCorner(bottomLeftCorner),
        _topRightCorner(autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength)),
        _density(density),
        _structure(structure),
        _centered(centered) {
    if (_particleSpacing <= 0.0 and _density > 0.0) {
      _particleSpacing = std::cbrt(std::sqrt(2.0) / _density);
    } else if (_particleSpacing > 0.0 and _density <= 0.0) {
      _density = std::sqrt(2.0) / (_particleSpacing * _particleSpacing * _particleSpacing);
    }
    _xOffset = _particleSpacing * 1. / 2.;
    _yOffset = _particleSpacing * sqrt(1. / 12.);
  }

  /**
   * Constructor with structure specified without density.
   * @param velocity
   * @param typeId
   * @param particleSpacing
   * @param boxLength
   * @param bottomLeftCorner
   * @param structure
   */
  CubeClosestPacked(const std::array<double, 3> &velocity, unsigned long typeId, double particleSpacing,
                    const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner,
                    Structure structure)
      : CubeClosestPacked(velocity, typeId, particleSpacing, boxLength, bottomLeftCorner, 0.0, structure) {}

  [[nodiscard]] double getParticleSpacing() const override { return _particleSpacing; }

  [[nodiscard]] double getParticleDensity() const { return _density; }

  [[nodiscard]] Structure getStructure() const { return _structure; }

  /**
   * Returns the total number of particles which will be / have been generated.
   * @return number of generated particles.
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    if (_structure == Structure::fcc) {
      return autopasTools::generators::FCCGenerator::getNumberOfParticles(_bottomLeftCorner, _topRightCorner,
                                                                          _particleSpacing);
    }

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
    const double numLayers = std::ceil(_boxLength[2] / spacingLayer);
    // Add up all even and odd layers.
    return evenLayer * std::ceil(numLayers / 2.) + oddLayer * std::floor(numLayers / 2.);
  }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    using namespace autopas::utils::ArrayMath::literals;
    return _bottomLeftCorner + _boxLength;
  }

  /**
   * Converts the object to a human-readable string
   * @return the generated string
   */
  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << _particleSpacing << "\n";
    output << std::setw(_valueOffset) << std::left << "particle-density"
           << ":  " << _density << "\n";
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << "\n";
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << "\n";
    output << std::setw(_valueOffset) << std::left << "structure"
           << ":  " << (_structure == Structure::fcc ? "fcc" : "hcp") << "\n";
    output << std::setw(_valueOffset) << std::left << "centered"
           << ":  " << std::to_string(_centered) << "\n";
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

    if (_structure == Structure::fcc) {
      autopasTools::generators::FCCGenerator::fillWithParticles(particlesWrapper, _bottomLeftCorner, _topRightCorner,
                                                                dummyParticle, _particleSpacing);
    } else {
      autopasTools::generators::HCPGenerator::fillWithParticles(particlesWrapper, _bottomLeftCorner, _topRightCorner,
                                                                dummyParticle, _particleSpacing);
    }
  }

 private:
  /**
   * The distance between the particles.
   */
  double _particleSpacing{1.0};

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
   * Shorter part of the bisectrix when split at the intersection of all bisectrices.
   */
  double _xOffset{0.5};

  /**
   * Shorter part of the bisectrix when split at the intersection of all bisectrices.
   */
  double _yOffset{0.288675};

  /**
   * Target particle density.
   */
  double _density{0.0};

  /**
   * Structure (fcc or hcp).
   */
  Structure _structure{Structure::fcc};

  /**
   * Lattice alignment (First particle at center or origin of a lattice unit cell).
   */
  bool _centered{true};
};
