/**
 * @file GaussianGenerator.h
 * @author F. Gratl
 * @date 7/31/18
 */

#pragma once

#include <random>
#include "autopas/AutoPas.h"

class GaussianGenerator {
 public:
  /**
   * fills a autopas object with gaussian distributed particles
   * @tparam Particle Type of particle to be generated
   * @tparam ParticleCell
   * @param autoPas
   * @param numParticles number of particles
   * @param defaultParticle inserted particle
   * @param distributionMean mean value / expected value
   * @param distributionStdDev standard deviation
   */
  template <class Particle, class ParticleCell>
  static void fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, size_t numParticles,
                                const Particle &defaultParticle = autopas::Particle(), double distributionMean = 5.0,
                                double distributionStdDev = 2.0);
};

template <class Particle, class ParticleCell>
void GaussianGenerator::fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, size_t numParticles,
                                          const Particle &defaultParticle, double distributionMean,
                                          double distributionStdDev) {
  std::default_random_engine generator(42);
  std::normal_distribution<double> distribution(distributionMean, distributionStdDev);

  for (size_t id = 0; id < numParticles;) {
    std::array<double, 3> position = {distribution(generator), distribution(generator), distribution(generator)};
    // only increment loop var (and place particle) if position is valid
    if (not autopas::utils::inBox(position, autoPas.getBoxMin(), autoPas.getBoxMax())) continue;
    Particle p(defaultParticle);
    p.setR(position);
    p.setID(id++);
    autoPas.addParticle(p);
  }
}