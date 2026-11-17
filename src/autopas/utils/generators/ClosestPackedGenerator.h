/**
 * @file ClosestPackedGenerator.h
 * @author S. Newcome
 * @date 28.09.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ParticleTypeTrait.h"
#include "autopas/utils/generators/FCCGenerator.h"
#include "autopas/utils/generators/HCPGenerator.h"

/**
 * Generator for closest packed particle lattices (FCC or HCP), including the choice of the particle spacing for a
 * target density. Moved here from md-flexible's CubeClosestPacked.
 */
namespace autopas::generators::ClosestPackedGenerator {

/**
 * Structure type of closest packing.
 */
enum LatticeStructure { FCC, HCP };

/**
 * Helper to calculate the total particle count for a given lattice structure and spacing.
 * @param boxMin
 * @param boxMax
 * @param spacing
 * @param structure
 * @param centered
 * @return particle count
 */
inline size_t calculateParticleCount(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                     const double spacing, const LatticeStructure structure, const bool centered) {
  switch (structure) {
    case FCC:
      return autopas::generators::FCCGenerator::getNumberOfParticles(boxMin, boxMax, spacing, centered);
    case HCP:
      return autopas::generators::HCPGenerator::getNumberOfParticles(boxMin, boxMax, spacing, centered);
    default:
      autopas::utils::ExceptionHandler::exception(
          "ClosestPackedGenerator: Unknown lattice structure. Possible values: (fcc hcp)");
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
inline double optimizeSpacingForDensity(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
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
 * Fills any container with particles arranged in a closest packed lattice.
 * Particle IDs start from the default particle.
 * @tparam Container Arbitrary container class that supports addParticle().
 * @param container
 * @param boxMin
 * @param boxMax
 * @param defaultParticle Blueprint particle.
 * @param spacing closest distance between neighboring particles
 * @param structure
 * @param centered If true, the cell offset moved inward by 1/4 * lattice constant.
 */
template <class Container>
void fillWithParticles(Container &container, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                       const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle,
                       const double spacing, const LatticeStructure structure, const bool centered) {
  switch (structure) {
    case FCC:
      autopas::generators::FCCGenerator::fillWithParticles(container, boxMin, boxMax, defaultParticle, spacing,
                                                           centered);
      break;
    case HCP:
      autopas::generators::HCPGenerator::fillWithParticles(container, boxMin, boxMax, defaultParticle, spacing,
                                                           centered);
      break;
    default:
      autopas::utils::ExceptionHandler::exception(
          "ClosestPackedGenerator: Unknown lattice structure. Possible values: (fcc hcp)");
  }
}

}  // namespace autopas::generators::ClosestPackedGenerator