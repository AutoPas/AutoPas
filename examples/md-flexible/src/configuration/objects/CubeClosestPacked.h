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
                    const LatticeStructure structure = HCP, const bool centered = false)
      : CubeClosestPacked(velocity, typeId, boxLength, bottomLeftCorner, structure, centered) {
    _particleSpacing = particleSpacing;
    _density =
        static_cast<double>(CubeClosestPacked::getParticlesTotal()) / (boxLength[0] * boxLength[1] * boxLength[2]);
  }

  /**
   * Constructor based on a given density. Spacing is automatically optimized (s >= s0) so total density is as close
   * to target as possible without compressing particles.
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
                    const LatticeStructure structure = HCP, const bool centered = false)
      : CubeClosestPacked(velocity, typeId, boxLength, bottomLeftCorner, structure, centered) {
    _particleSpacing = optimizeSpacingForDensity(_bottomLeftCorner, _topRightCorner, density, structure, centered);
    const auto newDensity =
        static_cast<double>(CubeClosestPacked::getParticlesTotal()) / (boxLength[0] * boxLength[1] * boxLength[2]);
    _density = newDensity;
    if (std::abs(newDensity - density) > 1e-10) {
      std::cout << "CubeClosestPacked: The requested density of " << density
                << " could not be achieved with the given box length. Actual density: " << newDensity << "."
                << std::endl;
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
    return calculateParticleCount(_bottomLeftCorner, _topRightCorner, _particleSpacing, _structure, _centered);
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
                    const bool centered = false)
      : Object(velocity, typeId),
        _boxLength(boxLength),
        _bottomLeftCorner(bottomLeftCorner),
        _topRightCorner(autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength)),
        _structure(structure),
        _centered(centered) {}

  /**
   * Helper to calculate the total particle count for a given lattice structure and spacing.
   * @param boxMin
   * @param boxMax
   * @param spacing
   * @param structure
   * @param centered
   * @return particle count
   */
  static size_t calculateParticleCount(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                       const double spacing, const LatticeStructure structure, const bool centered) {
    switch (structure) {
      case FCC:
        return autopas::generators::FCCGenerator::getNumberOfParticles(boxMin, boxMax, spacing, centered);
      case HCP:
        return autopas::generators::HCPGenerator::getNumberOfParticles(boxMin, boxMax, spacing, centered);
      default:
        autopas::utils::ExceptionHandler::exception(
            "CubeClosestPacked: Unknown lattice structure. Possible values: (fcc hcp)");
    }
    return 0;
  }

  /**
   * Computes an optimized particle spacing for a target density in a given bounding box.
   * Ensures particles are never compressed below the theoretical bulk spacing s0 = cbrt(sqrt(2) / density),
   * but may increase the spacing (s >= s0) to make the resulting total particle count as close to
   * round(density * volume) as possible.
   *
   * @param boxMin Minimum box coordinates.
   * @param boxMax Maximum box coordinates.
   * @param targetDensity Target particle density.
   * @param structure Lattice structure (FCC or HCP).
   * @param centered If true, lattice is centered.
   * @return Optimized particle spacing s >= s0.
   */
  static double optimizeSpacingForDensity(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                          const double targetDensity, const LatticeStructure structure,
                                          const bool centered) {
    const double volume = (boxMax[0] - boxMin[0]) * (boxMax[1] - boxMin[1]) * (boxMax[2] - boxMin[2]);

    // Spacing for a theoretical infinite lattice (s0)
    const double spacingExact = std::cbrt(std::sqrt(2.0) / targetDensity);
    const auto targetParticles = static_cast<size_t>(std::round(targetDensity * volume));

    const auto countAtExactSpacing = calculateParticleCount(boxMin, boxMax, spacingExact, structure, centered);

    if (countAtExactSpacing <= targetParticles) {
      return spacingExact;
    }

    // If countAtExactSpacing > targetParticles. We search for a spacing s >= s0 to reduce particle count to
    // targetParticles.
    double spacingLow = spacingExact;
    double spacingHigh = spacingExact * 1.2;
    size_t countAtHighSpacing = countAtExactSpacing;
    while (countAtHighSpacing > targetParticles) {
      countAtHighSpacing = calculateParticleCount(boxMin, boxMax, spacingHigh, structure, centered);
      spacingLow = spacingHigh;
      spacingHigh *= 1.2;
      if (countAtHighSpacing <= 1) {
        break;
      }
    }

    // Binary search for the transition around targetParticles
    constexpr size_t maxIterations = 20;
    constexpr double tolerance = 1e-10;
    for (size_t iter = 0; iter < maxIterations and (spacingHigh - spacingLow) > tolerance; ++iter) {
      const double spacingMid = spacingLow + 0.5 * (spacingHigh - spacingLow);
      const size_t countMid = calculateParticleCount(boxMin, boxMax, spacingMid, structure, centered);
      if (countMid > targetParticles) {
        spacingLow = spacingMid;
      } else {
        spacingHigh = spacingMid;
      }
    }

    // spacingLow gives count > targetParticles, spacingHigh gives count <= targetParticles.
    const size_t countLow = calculateParticleCount(boxMin, boxMax, spacingLow, structure, centered);
    const size_t countHigh = calculateParticleCount(boxMin, boxMax, spacingHigh, structure, centered);

    const size_t diffLow = (countLow >= targetParticles) ? (countLow - targetParticles) : (targetParticles - countLow);
    const size_t diffHigh =
        (countHigh >= targetParticles) ? (countHigh - targetParticles) : (targetParticles - countHigh);

    // Pick whichever is closer to targetParticles
    if (diffHigh <= diffLow) {
      return spacingHigh;
    }
    return spacingLow;
  }

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
  bool _centered{false};
};
