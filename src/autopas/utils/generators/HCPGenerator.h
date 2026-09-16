/**
 * @file HCPGenerator.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <array>
#include <cmath>

#include "autopas/utils/ParticleTypeTrait.h"

/**
 * Generator for an hexagonally closest packed particle grid.
 */
namespace autopas::generators::HCPGenerator {

/**
 * Calculates the number of particles generated in an HCP lattice within [boxMin, boxMax) in O(1) time.
 * @param boxMin
 * @param boxMax
 * @param spacing Nearest-neighbor distance d.
 * @param centeredAlignment If true, the cell offset moved inward by 1/4 * lattice constant.
 * @return Number of particles.
 */
inline size_t getNumberOfParticles(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                   const double spacing = 1.0, const bool centeredAlignment = true) {
  if (spacing <= 0.0) {
    return 0;
  }

  const double spacingRow = spacing * sqrt(3. / 4.);
  const double spacingLayer = spacing * sqrt(2. / 3.);
  const double xOffset = spacing * 1. / 2.;
  const double yOffset = spacing * sqrt(1. / 12.);

  std::array<double, 3> offset = {0.0, 0.0, 0.0};
  if (centeredAlignment) {
    offset = {spacing / 4.0, yOffset, spacingLayer / 2.0};
  }

  auto countSteps = [](const double start, const double stop, const double step) -> size_t {
    if (stop <= start or step <= 0.0) {
      return 0;
    }
    const auto count = static_cast<size_t>(std::ceil((stop - start) / step));
    return count;
  };

  const size_t numLayersZ = countSteps(boxMin[2] + offset[2], boxMax[2], spacingLayer);
  if (numLayersZ == 0) {
    return 0;
  }
  const size_t numEvenLayers = (numLayersZ + 1) / 2;
  const size_t numOddLayers = numLayersZ / 2;

  const double startYEven = boxMin[1] + offset[1];
  const size_t nyEvenLayer = countSteps(startYEven, boxMax[1], spacingRow);

  const double startYOdd = boxMin[1] + yOffset + offset[1];
  const size_t nyOddLayer = countSteps(startYOdd, boxMax[1], spacingRow);

  const double startXEvenRow = boxMin[0] + offset[0];
  const size_t nxEvenRow = countSteps(startXEvenRow, boxMax[0], spacing);

  const double startXOddRow = boxMin[0] + xOffset + offset[0];
  const size_t nxOddRow = countSteps(startXOddRow, boxMax[0], spacing);

  const size_t nyEvenLayerEvenRows = (nyEvenLayer + 1) / 2;
  const size_t nyEvenLayerOddRows = nyEvenLayer / 2;
  const size_t particlesPerEvenLayer = nyEvenLayerEvenRows * nxEvenRow + nyEvenLayerOddRows * nxOddRow;

  const size_t nyOddLayerOddRows = (nyOddLayer + 1) / 2;
  const size_t nyOddLayerEvenRows = nyOddLayer / 2;
  const size_t particlesPerOddLayer = nyOddLayerOddRows * nxOddRow + nyOddLayerEvenRows * nxEvenRow;

  return numEvenLayers * particlesPerEvenLayer + numOddLayers * particlesPerOddLayer;
}

/**
 * Fills any container (also AutoPas object) with hexagonally closest packed particles.
 * Particle properties will be used from the default particle. Particle IDs start from the default particle.
 * @tparam Container Arbitrary container class that needs to support addParticle().
 * @param container
 * @param boxMin
 * @param boxMax
 * @param defaultParticle
 * @param spacing Distance between all neighboring particles
 * @param centeredAlignment If true, the cell offset moved inward by 1/4 * lattice constant.
 */
template <class Container>
void fillWithParticles(Container &container, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                       const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle =
                           typename autopas::utils::ParticleTypeTrait<Container>::value(),
                       const double spacing = 1., const bool centeredAlignment = true) {
  // Spacing in y direction when only moving 60° on the unit circle. Or the height in an equilateral triangle.
  const double spacingRow = spacing * sqrt(3. / 4.);
  // Spacing in z direction. Height in an equilateral tetrahedron.
  const double spacingLayer = spacing * sqrt(2. / 3.);
  // Shorter part of the bisectrix when split at the intersection of all bisectrices.
  const double xOffset = spacing * 1. / 2.;
  // Shorter part of the bisectrix when split at the intersection of all bisectrices.
  const double yOffset = spacing * sqrt(1. / 12.);

  // The packing alternates between odd and even layers and rows
  bool evenLayer = true;

  // centered offset
  std::array<double, 3> offset = {0.0, 0.0, 0.0};
  if (centeredAlignment) {
    offset = {spacing / 4.0, yOffset, spacingLayer / 2.0};
  }

  size_t id = defaultParticle.getID();
  for (double z = boxMin[2] + offset[2]; z < boxMax[2]; z += spacingLayer) {
    const double startY = (evenLayer ? boxMin[1] : boxMin[1] + yOffset) + offset[1];
    bool evenRow = evenLayer;  // To ensure layers are alternating as for hexagonal close packed.
    for (double y = startY; y < boxMax[1]; y += spacingRow) {
      const double startX = (evenRow ? boxMin[0] : boxMin[0] + xOffset) + offset[0];
      for (double x = startX; x < boxMax[0]; x += spacing) {
        auto p = defaultParticle;
        p.setR({x, y, z});
        p.setID(id++);
        container.addParticle(p);
      }
      evenRow = not evenRow;
    }
    evenLayer = not evenLayer;
  }
};
};  // namespace autopas::generators::HCPGenerator
