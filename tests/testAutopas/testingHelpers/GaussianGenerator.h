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
  template <class Particle, class ParticleCell, typename precision = double>
  static void fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, size_t numParticles,
                                const Particle &defaultParticle = autopas::Particle(), precision distributionMean = 5.0,
                                precision distributionStdDev = 2.0);
};

template <class Particle, class ParticleCell, typename precision = double>
void GaussianGenerator::fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, size_t numParticles,
                                          const Particle &defaultParticle, precision distributionMean,
                                          precision distributionStdDev) {
  std::default_random_engine generator(42);
  std::normal_distribution<precision> distribution(distributionMean, distributionStdDev);

  for (size_t id = 0; id < numParticles;) {
    std::array<precision, 3> position = {distribution(generator), distribution(generator), distribution(generator)};
    // only increment loop var (and place particle) if position is valid
    if (not autopas::utils::inBox(position, autoPas.getBoxMin(), autoPas.getBoxMax())) continue;
    Particle p(defaultParticle);
    p.setR(position);
    p.setID(id++);
    autoPas.addParticle(p);
  }
}
