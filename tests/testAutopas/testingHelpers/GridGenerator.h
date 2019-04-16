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
   * fills a cell vector with a cuboid mesh of particles
   * @tparam Particle Type of particle to be generated
   * @tparam ParticleCell
   * @param cells
   * @param particlesPerDim Number of particles per dimension
   * @param defaultParticle inserted particle
   * @param spacing Factor for distance between two particles along one
   * dimension (default is 1)
   * @param offset Offset to move all particles
   */
  template <class Particle, class ParticleCell>
  static void fillWithParticles(std::vector<ParticleCell>& cells, const std::array<size_t, 3>& particlesPerDim,
                                const Particle& defaultParticle = autopas::Particle(),
                                const std::array<double, 3>& spacing = {1.0, 1.0, 1.0},
                                const std::array<double, 3>& offset = {.5, .5, .5});

  /**
   * fills a autopas object with a cuboid mesh of particles
   * @tparam Particle Type of particle to be generated
   * @tparam ParticleCell
   * @param autoPas
   * @param particlesPerDim Number of particles per dimension.
   * @param defaultParticle inserted particle
   * @param spacing Factor for distance between two particles along one
   * dimension (default is 1).
   * @param offset Offset to move all particles.
   */
  template <class Particle, class ParticleCell>
  static void fillWithParticles(autopas::AutoPas<Particle, ParticleCell>& autoPas,
                                const std::array<size_t, 3>& particlesPerDim,
                                const Particle& defaultParticle = autopas::Particle(),
                                const std::array<double, 3>& spacing = {1.0, 1.0, 1.0},
                                const std::array<double, 3>& offset = {.5, .5, .5});
};

template <class Particle, class ParticleCell>
void GridGenerator::fillWithParticles(std::vector<ParticleCell>& cells, const std::array<size_t, 3>& particlesPerDim,
                                      const Particle& defaultParticle, const std::array<double, 3>& spacing,
                                      const std::array<double, 3>& offset) {
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

template <class Particle, class ParticleCell>
void GridGenerator::fillWithParticles(autopas::AutoPas<Particle, ParticleCell>& autoPas,
                                      const std::array<size_t, 3>& particlesPerDim, const Particle& defaultParticle,
                                      const std::array<double, 3>& spacing, const std::array<double, 3>& offset) {
  size_t id = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        Particle p(defaultParticle);
        p.setR({x * spacing[0] + offset[0], y * spacing[1] + offset[1], z * spacing[2] + offset[2]});
        p.setID(id++);
        autoPas.addParticle(p);
      }
    }
  }
}