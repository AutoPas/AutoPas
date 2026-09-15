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
namespace autopasTools::generators::HCPGenerator {

/**
 * Calculates the number of particles generated in an HCP lattice within [boxMin, boxMax).
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

  size_t count = 0;
  bool evenLayer = true;
  for (double z = boxMin[2] + offset[2]; z < boxMax[2]; z += spacingLayer) {
    const double startY = (evenLayer ? boxMin[1] : boxMin[1] + yOffset) + offset[1];
    bool evenRow = evenLayer;
    for (double y = startY; y < boxMax[1]; y += spacingRow) {
      const double startX = (evenRow ? boxMin[0] : boxMin[0] + xOffset) + offset[0];
      for (double x = startX; x < boxMax[0]; x += spacing) {
        ++count;
      }
      evenRow = not evenRow;
    }
    evenLayer = not evenLayer;
  }
  return count;
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
};  // namespace autopasTools::generators::HCPGenerator
