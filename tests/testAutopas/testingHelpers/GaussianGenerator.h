/**
 * @file GaussianGenerator.h
 * @author F. Gratl
 * @date 7/31/18
 */

#pragma once

#include <random>
#include "autopas/AutoPas.h"

/**
 * Generator class for gaussian distributions
 */
class GaussianGenerator {
 public:
  /**
   * Fills any container (also AutoPas object) with randomly 3D gaussian distributed particles.
   *
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param autoPas
   * @param numParticles number of particles
   * @param defaultParticle inserted particle
   * @param distributionMean mean value / expected value
   * @param distributionStdDev standard deviation
   */
  template <class Particle, class ParticleCell>
  static void fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, size_t numParticles,
                                const Particle &defaultParticle = autopas::Particle(),
                                typename Particle::ParticleFloatingPointType distributionMean = 5.0,
                                typename Particle::ParticleFloatingPointType distributionStdDev = 2.0);
};

template <class Particle, class ParticleCell>
void GaussianGenerator::fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, size_t numParticles,
                                          const Particle &defaultParticle,
                                          typename Particle::ParticleFloatingPointType distributionMean,
                                          typename Particle::ParticleFloatingPointType distributionStdDev) {
  std::default_random_engine generator(42);
  std::normal_distribution<typename Particle::ParticleFloatingPointType> distribution(distributionMean,
                                                                                      distributionStdDev);

  for (size_t id = 0; id < numParticles;) {
    std::array<typename Particle::ParticleFloatingPointType, 3> position = {
        distribution(generator), distribution(generator), distribution(generator)};
    // only increment loop var (and place particle) if position is valid
    if (not autopas::utils::inBox(position, autoPas.getBoxMin(), autoPas.getBoxMax())) continue;
    Particle p(defaultParticle);
    p.setR(position);
    p.setID(id++);
    autoPas.addParticle(p);
  }
}
