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
   * Particle properties will be used from the default particle. Particle IDs start from the default particle.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param autoPas
   * @param boxMin
   * @param boxMax
   * @param numParticles number of particles
   * @param defaultParticle inserted particle
   * @param distributionMean mean value / expected value
   * @param distributionStdDev standard deviation
   */
  template <class Particle, class ParticleCell>
  static void fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, const std::array<double, 3> &boxMin,
                                const std::array<double, 3> &boxMax, size_t numParticles,
                                const Particle &defaultParticle = autopas::MoleculeLJ<>(),
                                double distributionMean = 5.0, double distributionStdDev = 2.0);
};

template <class Particle, class ParticleCell>
void GaussianGenerator::fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas,
                                          const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                          size_t numParticles, const Particle &defaultParticle, double distributionMean,
                                          double distributionStdDev) {
  std::default_random_engine generator(42);
  std::normal_distribution<double> distribution(distributionMean, distributionStdDev);

  for (unsigned long i = defaultParticle.getID(); i < defaultParticle.getID() + numParticles; ++i) {
    std::array<double, 3> position;
    // assert that position is valid
    do {
      position = {distribution(generator), distribution(generator), distribution(generator)};
    } while (not autopas::utils::inBox(position, boxMin, boxMax));
    Particle p(defaultParticle);
    p.setR(position);
    p.setID(i);
    autoPas.addParticle(p);
  }
}