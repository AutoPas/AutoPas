/**
 * @file ClosestPackingGenerator.h
 * @author F. Gratl
 * @date 21.08.20
 */

#pragma once

#include <cmath>
#include <vector>

#include "autopas/utils/ParticleTypeTrait.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

/**
 * Generator for grids of particles.
 */
namespace autopasTools::generators::ClosestPackingGenerator {
#if AUTOPAS_PRECISION_MODE == DPDP
using CalcPrecision = double;
#else
using CalcPrecision = float;
#endif
/**
 * Fills any container (also AutoPas object) with a hexagonally closest packed particles.
 * Particle properties will be used from the default particle. Particle IDs start from the default particle.
 * @tparam Container Arbitrary container class that needs to support addParticle().
 * @param container
 * @param boxMin
 * @param boxMax
 * @param defaultParticle
 * @param spacing Distance between all neighboring particles
 */
template <class Container>
void fillWithParticles(Container &container, const std::array<CalcPrecision, 3> &boxMin,
                       const std::array<CalcPrecision, 3> &boxMax,
                       const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle =
                           typename autopas::utils::ParticleTypeTrait<Container>::value(),
                       const CalcPrecision spacing = 1.) {
  // Spacing in y direction when only moving 60° on the unit circle. Or the height in an equilateral triangle.
  const CalcPrecision spacingRow = spacing * sqrt(3. / 4.);
  // Spacing in z direction. Height in an equilateral tetraeder.
  const CalcPrecision spacingLayer = spacing * sqrt(2. / 3.);
  // Shorter part of the bisectrix when split at the intersection of all bisectrices.
  const CalcPrecision xOffset = spacing * 1. / 2.;
  // Shorter part of the bisectrix when split at the intersection of all bisectrices.
  const CalcPrecision yOffset = spacing * sqrt(1. / 12.);

  // The packing alternates between odd and even layers and rows
  bool evenLayer = true;

  size_t id = defaultParticle.getID();
  for (CalcPrecision z = boxMin[2]; z < boxMax[2]; z += spacingLayer) {
    CalcPrecision starty = evenLayer ? boxMin[1] : boxMin[1] + yOffset;
    bool evenRow = evenLayer;  // To ensure layers are alternating as for hexagonal close packed.
    for (CalcPrecision y = starty; y < boxMax[1]; y += spacingRow) {
      CalcPrecision startx = evenRow ? boxMin[0] : boxMin[0] + xOffset;
      for (CalcPrecision x = startx; x < boxMax[0]; x += spacing) {
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
};  // namespace autopasTools::generators::ClosestPackingGenerator
