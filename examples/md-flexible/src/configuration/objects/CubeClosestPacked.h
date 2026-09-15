/**
 * @file CubeClosestPacked.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <cmath>

#include "Object.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/generators/FCCGenerator.h"
#include "autopas/utils/generators/HCPGenerator.h"
#include "autopas/utils/generators/PseudoContainer.h"

/**
 * Class describing a cube of closest packed particles (FCC or HCP).
 */
class CubeClosestPacked : public Object {
 public:
  /**
   * Structure type of closest packing.
   */
  enum LatticeStructure { FCC, HCP };

  /**
   * Constructor based on a given particle spacing.
   * @param velocity
   * @param typeId
   * @param particleSpacing closest distance between neighboring particles
   * @param boxLength
   * @param bottomLeftCorner
   * @param structure fcc (default) or hcp.
   * @param centered If true, the cell offset moved inward by 1/4 * lattice constant.
   */
  CubeClosestPacked(const std::array<double, 3> &velocity, const size_t typeId, const double particleSpacing,
                    const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner,
                    const LatticeStructure structure = HCP, const bool centered = true)
      : CubeClosestPacked(velocity, typeId, boxLength, bottomLeftCorner, structure, centered) {
    _particleSpacing = particleSpacing;
    _density =
        static_cast<double>(CubeClosestPacked::getParticlesTotal()) / (boxLength[0] * boxLength[1] * boxLength[2]);
  }

  /**
   * Constructor based on a given density.
   * @param velocity
   * @param typeId
   * @param boxLength
   * @param bottomLeftCorner
   * @param density Target particle density.
   * @param structure hcp (default) or fcc.
   * @param centered If true, the cell offset moved inward by 1/4 * lattice constant.
   */
  CubeClosestPacked(const std::array<double, 3> &velocity, const size_t typeId, const std::array<double, 3> &boxLength,
                    const std::array<double, 3> &bottomLeftCorner, const double density,
                    const LatticeStructure structure = HCP, const bool centered = true)
      : CubeClosestPacked(velocity, typeId, boxLength, bottomLeftCorner, structure, centered) {
    _particleSpacing = std::cbrt(std::sqrt(2.0) / density);
    const auto newDensity =
        static_cast<double>(CubeClosestPacked::getParticlesTotal()) / (boxLength[0] * boxLength[1] * boxLength[2]);
    _density = newDensity;
    if (std::abs(newDensity - density) > 1e-10) {
      std::cout << "CubeClosestPacked: The requested density of " << density
                << " could not be achieved with the given box "
                << "length. Actual density: " << newDensity << "." << std::endl;
    }
  }

  [[nodiscard]] std::string getObjectType() const override { return "CubeClosestPacked"; }

  [[nodiscard]] double getParticleSpacing() const override { return _particleSpacing; }

  /**
   * Returns the particle density.
   * @return density of particles.
   */
  [[nodiscard]] double getParticleDensity() const { return _density; }

  /**
   * Returns the lattice structure of the closest packing generator.
   * @return the lattice structure of the closest packing generator.
   */
  [[nodiscard]] LatticeStructure getLatticeStructure() const { return _structure; }

  /**
   * Returns the total number of particles which will be / have been generated.
   * @return number of generated particles.
   */
  [[nodiscard]] size_t getParticlesTotal() const override {
    switch (_structure) {
      case FCC:
        return autopasTools::generators::FCCGenerator::getNumberOfParticles(_bottomLeftCorner, _topRightCorner,
                                                                            _particleSpacing, _centered);
      case HCP:
        return autopasTools::generators::HCPGenerator::getNumberOfParticles(_bottomLeftCorner, _topRightCorner,
                                                                            _particleSpacing, _centered);
      default:
        autopas::utils::ExceptionHandler::exception(
            "CubeClosestPacked: Unknown lattice structure. Possible values: (fcc hcp)");
    }
    return 0;
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
           << ":  " << (_structure == FCC ? "fcc" : "hcp") << "\n";
    output << std::setw(_valueOffset) << std::left << "centered"
           << ":  " << (_centered ? "true" : "false") << "\n";
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates particles based on the parameters provided to the CubeClosestPacked Object in the configuration file.
   * @param particles: The container, where the new particles get stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    // Wrapper so that std::vector can be used as an AutoPas::ParticleContainer
    auto particlesWrapper = autopas::generators::PseudoContainer(particles);

    // dummy particle used as a template with id of the first newly generated one
    const ParticleType dummyParticle = getDummyParticle(particles.size());

    switch (_structure) {
      case FCC:
        autopas::generators::FCCGenerator::fillWithParticles(particlesWrapper, _bottomLeftCorner, _topRightCorner,
                                                                  dummyParticle, _particleSpacing, _centered);
        break;
      case HCP:
        autopas::generators::HCPGenerator::fillWithParticles(particlesWrapper, _bottomLeftCorner, _topRightCorner,
                                                                  dummyParticle, _particleSpacing, _centered);
        break;
      default:
        autopas::utils::ExceptionHandler::exception(
            "CubeClosestPacked: Unknown lattice structure. Possible values: (fcc hcp)");
    }
  }

 private:
  /**
   * Internal constructor.
   * @param velocity
   * @param typeId
   * @param boxLength
   * @param bottomLeftCorner
   * @param structure fcc (default) or hcp.
   * @param centered If true, the cell offset moved inward by 1/4 * lattice constant.
   */
  CubeClosestPacked(const std::array<double, 3> &velocity, const size_t typeId, const std::array<double, 3> &boxLength,
                    const std::array<double, 3> &bottomLeftCorner, const LatticeStructure structure = FCC,
                    const bool centered = true)
      : Object(velocity, typeId),
        _boxLength(boxLength),
        _bottomLeftCorner(bottomLeftCorner),
        _topRightCorner(autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength)),
        _structure(structure),
        _centered(centered) {}

  /**
   * The distance between the particles.
   */
  double _particleSpacing{0.0};

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
   * Target particle density.
   */
  double _density{0.0};

  /**
   * Structure (hcp or fcc).
   */
  LatticeStructure _structure{HCP};

  /**
   * Lattice alignment (First particle at center or origin of a lattice unit cell).
   */
  bool _centered{true};
};
