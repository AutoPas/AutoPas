/**
 * @file GaussianGenerator.h
 * @author F. Gratl
 * @date 5/25/18
 */

#pragma once

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ParticleTypeTrait.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopasTools::generators {
/**
 * Generator for grids of particles.
 */
namespace GridGenerator {
/**
 * Fills a cell vector with a cuboid mesh of particles.
 *
 * @tparam ParticleCell Type of ParticleCells in cell vector.
 * @param cells Cell vector.
 * @param cellsPerDimension Number of cells per dimension.
 * @param particlesPerDim Number of particles per dimension.
 * @param defaultParticle
 * @param spacing Factor for distance between two particles along one dimension (default is 1).
 * @param offset Offset to move all particles.
 * @param cellSize The size of a cell.
 */
template <class ParticleCell>
void fillWithParticles(
    std::vector<ParticleCell> &cells, const std::array<size_t, 3> &cellsPerDimension,
    const std::array<size_t, 3> &particlesPerDim,
    const typename ParticleCell::ParticleType &defaultParticle = typename ParticleCell::ParticleType(),
    const std::array<CalcPrecision, 3> &spacing = std::array<CalcPrecision, 3>{1, 1, 1},
    const std::array<CalcPrecision, 3> &offset = std::array<CalcPrecision, 3>{.5, .5, .5},
    const std::array<CalcPrecision, 3> &cellSize = {1., 1., 1.});

/**
 * Fills any container (also AutoPas object) with a cuboid mesh of particles.
 * Particle properties will be used from the default particle. Particle IDs start from the default particle.
 * @tparam Container Arbitrary container class that needs to support addParticle().
 * @param container
 * @param particlesPerDim Number of particles per dimension.
 * @param defaultParticle
 * @param spacing Factor for distance between two particles along one dimension (default is 1).
 * @param offset Offset to move all particles.
 */
template <class Container>
void fillWithParticles(Container &container, const std::array<size_t, 3> &particlesPerDim,
                       const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle =
                           typename autopas::utils::ParticleTypeTrait<Container>::value(),
                       const std::array<CalcPrecision, 3> &spacing = std::array<CalcPrecision, 3>{1., 1., 1.},
                       const std::array<CalcPrecision, 3> &offset = std::array<CalcPrecision, 3>{.5, .5, .5});
};  // namespace GridGenerator

template <class ParticleCell>
void GridGenerator::fillWithParticles(std::vector<ParticleCell> &cells, const std::array<size_t, 3> &cellsPerDimension,
                                      const std::array<size_t, 3> &particlesPerDim,
                                      const typename ParticleCell::ParticleType &defaultParticle,
                                      const std::array<CalcPrecision, 3> &spacing,
                                      const std::array<CalcPrecision, 3> &offset,
                                      const std::array<CalcPrecision, 3> &cellSize) {
  size_t id = defaultParticle.getID();
  for (unsigned long z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned long y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned long x = 0; x < particlesPerDim[0]; ++x) {
        auto p = defaultParticle;
        std::array<CalcPrecision, 3> pos{static_cast<CalcPrecision>(x) * spacing[0] + offset[0],
                                         static_cast<CalcPrecision>(y) * spacing[1] + offset[1],
                                         static_cast<CalcPrecision>(z) * spacing[2] + offset[2]};
        std::array<unsigned long, 3> cellIndex3D{static_cast<unsigned long>(pos[0] / cellSize[0]),
                                                 static_cast<unsigned long>(pos[1] / cellSize[1]),
                                                 static_cast<unsigned long>(pos[2] / cellSize[2])};
        p.setR(pos);
        p.setID(id++);
        const auto cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(cellIndex3D, cellsPerDimension);

        if (cells[cellIndex].getPossibleParticleOwnerships() == autopas::OwnershipState::halo) {
          p.setOwnershipState(autopas::OwnershipState::halo);
        }
        if (cells[cellIndex].getPossibleParticleOwnerships() == autopas::OwnershipState::dummy) {
          p.setOwnershipState(autopas::OwnershipState::dummy);
        }
        cells[cellIndex].addParticle(p);
      }
    }
  }
}

template <class Container>
void GridGenerator::fillWithParticles(
    Container &container, const std::array<size_t, 3> &particlesPerDim,
    const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle,
    const std::array<CalcPrecision, 3> &spacing, const std::array<CalcPrecision, 3> &offset) {
  size_t id = defaultParticle.getID();
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        auto p = defaultParticle;
        p.setR({x * spacing[0] + offset[0], y * spacing[1] + offset[1], z * spacing[2] + offset[2]});
        p.setID(id++);
        container.addParticle(p);
      }
    }
  }
}
}  // namespace autopasTools::generators
