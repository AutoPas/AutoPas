/**
 * GridGenerator.h
 *
 *  Created on: 5/25/18
 *     Aauthor: F. Gratl
 */

#pragma once
#include "AutoPas.h"

/**
 * Generator for grids of particles.
 */
class GridGenerator {
 public:

/**
 * fills a cell vector with a cuboid mesh of particles
 * @tparam Particle Type of particle to be generated
 * @tparam ParticleCell
 * @param cells
 * @param particlesPerDim Number of particles per dimension
 * @param spacing Factor for distance between two particles along one dimension (default is 1)
 * @param offset Offset to move all particles
 */
  template<class Particle, class ParticleCell>
  static void fillWithParticles(std::vector<ParticleCell> &cells,
                                std::array<size_t, 3> particlesPerDim,
                                Particle &defaultParicle = autopas::Particle(),
                                std::array<double, 3> spacing = std::array<double, 3>{1, 1, 1},
                                std::array<double, 3> offset = std::array<double, 3>{.5, .5, .5});

/**
 * fills a autopas object with a cuboid mesh of particles
 * @tparam Particle Type of particle to be generated
 * @tparam ParticleCell
 * @param autoPas
 * @param particlesPerDim Number of particles per dimension.
 * @param spacing Factor for distance between two particles along one dimension (default is 1).
 * @param offset Offset to move all particles.
 */
  template<class Particle, class ParticleCell>
  static void fillWithParticles(AutoPas<Particle, ParticleCell> &autoPas,
                                std::array<size_t, 3> particlesPerDim,
                                Particle &defaultParicle = autopas::Particle(),
                                std::array<double, 3> spacing = std::array<double, 3>{1, 1, 1},
                                std::array<double, 3> offset = std::array<double, 3>{.5, .5, .5});
};

template<class Particle, class ParticleCell>
void GridGenerator::fillWithParticles(
    std::vector<ParticleCell> &cells,
    std::array<size_t, 3> particlesPerDim,
    Particle &defaultParicle,
    std::array<double, 3> spacing,
    std::array<double, 3> offset) {
  size_t id = 0;
  size_t cellId = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        Particle p(defaultParicle);
        p.setR({x * spacing[0] + offset[0],
                y * spacing[1] + offset[1],
                z * spacing[2] + offset[2]});
        p.setID(id++);
        cells[cellId++].addParticle(p);
      }
    }
  }
}

template<class Particle, class ParticleCell>
void GridGenerator::fillWithParticles(AutoPas<Particle, ParticleCell> &autoPas,
                                                              std::array<size_t, 3> particlesPerDim,
                                                              Particle &defaultParicle,
                                                              std::array<double, 3> spacing,
                                                              std::array<double, 3> offset) {
  size_t id = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        Particle p(defaultParicle);
        p.setR({x * spacing[0] + offset[0],
                y * spacing[1] + offset[1],
                z * spacing[2] + offset[2]});
        p.setID(id++);
        autoPas.addParticle(p);
      }
    }
  }
}