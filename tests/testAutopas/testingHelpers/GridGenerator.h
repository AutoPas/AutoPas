/**
 * GridGenerator.h
 *
 *  Created on: 5/25/18
 *     Aauthor: F. Gratl
 */

#pragma once
#include "AutoPas.h"

template<class Particle, class ParticleCell>
class GridGenerator {
 public:
/**
 * fills a cell vector with a cuboid mesh of particles
 * @param cells
 * @param particlesPerDim number of particles per dimension
 * @param spacing factor for distance between two particles along one dimension (default is 1)
 * @param offset offset to move all particles
 */
  static void fillWithParticles(std::vector<ParticleCell> &cells,
                                std::array<size_t, 3> particlesPerDim,
                                std::array<double, 3> spacing = std::array<double, 3>{1, 1, 1},
                                std::array<double, 3> offset = std::array<double, 3>{.5, .5, .5});

/**
 * fills a autopas object with a cuboid mesh of particles
 * @param autoPas
 * @param particlesPerDim number of particles per dimension
 * @param spacing factor for distance between two particles along one dimension (default is 1)
 * @param offset offset to move all particles
 */
  static void fillWithParticles(AutoPas<Particle, ParticleCell> &autoPas,
                                std::array<size_t, 3> particlesPerDim,
                                std::array<double, 3> spacing = std::array<double, 3>{1, 1, 1},
                                std::array<double, 3> offset = std::array<double, 3>{.5, .5, .5});
};

template<class Particle, class ParticleCell>
void GridGenerator<Particle, ParticleCell>::fillWithParticles(
    std::vector<ParticleCell> &cells,
    std::array<size_t, 3> particlesPerDim,
    std::array<double, 3> spacing,
    std::array<double, 3> offset) {
  size_t id = 0;
  size_t cellId = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        auto p = autopas::Particle({x * spacing[0] + offset[0],
                                    y * spacing[1] + offset[1],
                                    z * spacing[2] + offset[2]},
                                   {0, 0, 0},
                                   id++);
        cells[cellId++].addParticle(p);
      }
    }
  }
}

template<class Particle, class ParticleCell>
void GridGenerator<Particle, ParticleCell>::fillWithParticles(AutoPas<Particle, ParticleCell> &autoPas,
                                                              std::array<size_t, 3> particlesPerDim,
                                                              std::array<double, 3> spacing,
                                                              std::array<double, 3> offset) {
  size_t id = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        auto p = autopas::MoleculeLJ({x * spacing[0] + offset[0],
                                      y * spacing[1] + offset[1],
                                      z * spacing[2] + offset[2]},
                                     {0, 0, 0}, id++);
        autoPas.addParticle(p);
      }
    }
  }
}