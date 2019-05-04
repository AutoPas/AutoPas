/**
 * @file GaussianGenerator.h
 * @author F. Gratl
 * @date 5/25/18
 */

#pragma once

#include "autopas/AutoPas.h"

/**
 * Generator for grids of particles.
 */
class GridGenerator {
 public:
  /**
   * Fills a cell vector with a cuboid mesh of particles.
   *
   * @tparam Particle Type of particle to be generated
   * @tparam ParticleCell
   * @param cells
   * @param particlesPerDim Number of particles per dimension
   * @param defaultParticle
   * @param spacing Factor for distance between two particles along one dimension (default is 1)
   * @param offset Offset to move all particles
   */
  template <class Particle, class ParticleCell>
  static void fillWithParticles(std::vector<ParticleCell> &cells, std::array<size_t, 3> particlesPerDim,
                                const Particle &defaultParticle = autopas::Particle(),
                                const std::array<double, 3> &spacing = std::array<double, 3>{1, 1, 1},
                                const std::array<double, 3> &offset = std::array<double, 3>{.5, .5, .5});

  /**
   * Fills any container (also AutoPas object) with randomly 3D gaussian distributed particles.
   *
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param particlesPerDim Number of particles per dimension.
   * @param defaultParticle
   * @param spacing Factor for distance between two particles along one dimension (default is 1).
   * @param offset Offset to move all particles.
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container &container, std::array<size_t, 3> particlesPerDim,
                                const Particle &defaultParticle = autopas::Particle(),
                                const std::array<double, 3> &spacing = std::array<double, 3>{1, 1, 1},
                                const std::array<double, 3> &offset = std::array<double, 3>{.5, .5, .5});
};

template <class Particle, class ParticleCell>
void GridGenerator::fillWithParticles(std::vector<ParticleCell> &cells, std::array<size_t, 3> particlesPerDim,
                                      const Particle &defaultParticle, const std::array<double, 3> &spacing,
                                      const std::array<double, 3> &offset) {
  size_t id = 0;
  size_t cellId = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        Particle p(defaultParticle);
        p.setR({x * spacing[0] + offset[0], y * spacing[1] + offset[1], z * spacing[2] + offset[2]});
        p.setID(id++);
        cells[cellId++].addParticle(p);
      }
    }
  }
}

template <class Container, class Particle>
void GridGenerator::fillWithParticles(Container &container, std::array<size_t, 3> particlesPerDim,
                                      const Particle &defaultParticle, const std::array<double, 3> &spacing,
                                      const std::array<double, 3> &offset) {
  size_t id = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        Particle p(defaultParticle);
        p.setR({x * spacing[0] + offset[0], y * spacing[1] + offset[1], z * spacing[2] + offset[2]});
        p.setID(id++);
        container.addParticle(p);
      }
    }
  }
}
