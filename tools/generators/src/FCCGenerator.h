/**
 * @file FCCGenerator.h
 * @author muehlhaeusser
 * @date 05.09.2026
 */

#pragma once

#include <array>
#include <cmath>

#include "autopas/utils/ParticleTypeTrait.h"

/**
 * Generator for an face-centered cubic (FCC) particle grid.
 */
namespace autopasTools::generators::FCCGenerator {

/**
 * Basis coordinates for conventional cubic FCC unit cell normalized by lattice constant a.
 */
inline constexpr std::array<std::array<double, 3>, 4> fccBasis = {{
    {0.0, 0.0, 0.0},
    {0.5, 0.5, 0.0},
    {0.5, 0.0, 0.5},
    {0.0, 0.5, 0.5},
}};

/**
 * Calculates the number of particles generated in an FCC lattice within [boxMin, boxMax).
 * @param boxMin
 * @param boxMax
 * @param spacing Nearest-neighbor distance d.
 * @param centeredAlignment If true, the cell offset moved inward by 1/4 * lattice constant.
 * @return Number of particles.
 */
inline size_t getNumberOfParticles(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                   const double spacing = 1.0, const bool centeredAlignment = true) {
  const double latticeConstant = std::sqrt(2.0) * spacing;
  if (spacing <= 0.0) {
    return 0;
  }

  std::array<double, 3> offset = {0.0, 0.0, 0.0};
  if (centeredAlignment) {
    offset = {latticeConstant / 4.0, latticeConstant / 4.0, latticeConstant / 4.0};
  }

  const size_t numCellsX = std::max<size_t>(1, std::ceil((boxMax[0] - boxMin[0]) / latticeConstant));
  const size_t numCellsY = std::max<size_t>(1, std::ceil((boxMax[1] - boxMin[1]) / latticeConstant));
  const size_t numCellsZ = std::max<size_t>(1, std::ceil((boxMax[2] - boxMin[2]) / latticeConstant));

  size_t count = 0;
  for (size_t z = 0; z < numCellsZ; ++z) {
    for (size_t y = 0; y < numCellsY; ++y) {
      for (size_t x = 0; x < numCellsX; ++x) {
        const std::array<double, 3> cellOrigin = {
            boxMin[0] + static_cast<double>(x) * latticeConstant + offset[0],
            boxMin[1] + static_cast<double>(y) * latticeConstant + offset[1],
            boxMin[2] + static_cast<double>(z) * latticeConstant + offset[2],
        };
        for (const auto &b : fccBasis) {
          const std::array<double, 3> pos = {
              cellOrigin[0] + b[0] * latticeConstant,
              cellOrigin[1] + b[1] * latticeConstant,
              cellOrigin[2] + b[2] * latticeConstant,
          };
          if (pos[0] >= boxMin[0] and pos[0] < boxMax[0] and pos[1] >= boxMin[1] and pos[1] < boxMax[1] and
              pos[2] >= boxMin[2] and pos[2] < boxMax[2]) {
            ++count;
          }
        }
      }
    }
  }
  return count;
}

/**
 * Fills any container with particles arranged in a conventional Face-Centered Cubic (FCC) lattice.
 * Particle IDs start from the default particle.
 * @tparam Container Arbitrary container class that supports addParticle().
 * @param container
 * @param boxMin
 * @param boxMax
 * @param defaultParticle Blueprint particle.
 * @param spacing Nearest-neighbor distance d.
 * @param centeredAlignment If true, the cell offset moved inward by 1/4 * lattice constant.
 */
template <class Container>
void fillWithParticles(Container &container, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                       const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle =
                           typename autopas::utils::ParticleTypeTrait<Container>::value(),
                       const double spacing = 1.0, const bool centeredAlignment = true) {
  const double latticeConstant = std::sqrt(2.0) * spacing;

  std::array<double, 3> offset = {0.0, 0.0, 0.0};
  if (centeredAlignment) {
    offset = {latticeConstant / 4.0, latticeConstant / 4.0, latticeConstant / 4.0};
  }

  const size_t numCellsX = std::max<size_t>(1, std::ceil((boxMax[0] - boxMin[0]) / latticeConstant));
  const size_t numCellsY = std::max<size_t>(1, std::ceil((boxMax[1] - boxMin[1]) / latticeConstant));
  const size_t numCellsZ = std::max<size_t>(1, std::ceil((boxMax[2] - boxMin[2]) / latticeConstant));

  size_t id = defaultParticle.getID();
  for (size_t z = 0; z < numCellsZ; ++z) {
    for (size_t y = 0; y < numCellsY; ++y) {
      for (size_t x = 0; x < numCellsX; ++x) {
        const std::array<double, 3> cellOrigin = {
            boxMin[0] + static_cast<double>(x) * latticeConstant + offset[0],
            boxMin[1] + static_cast<double>(y) * latticeConstant + offset[1],
            boxMin[2] + static_cast<double>(z) * latticeConstant + offset[2],
        };
        for (const auto &basis : fccBasis) {
          const std::array<double, 3> pos = {
              cellOrigin[0] + basis[0] * latticeConstant,
              cellOrigin[1] + basis[1] * latticeConstant,
              cellOrigin[2] + basis[2] * latticeConstant,
          };
          if (pos[0] >= boxMin[0] and pos[0] < boxMax[0] and pos[1] >= boxMin[1] and pos[1] < boxMax[1] and
              pos[2] >= boxMin[2] and pos[2] < boxMax[2]) {
            auto p = defaultParticle;
            p.setR(pos);
            p.setID(id++);
            container.addParticle(p);
          }
        }
      }
    }
  }
}

}  // namespace autopasTools::generators::FCCGenerator
