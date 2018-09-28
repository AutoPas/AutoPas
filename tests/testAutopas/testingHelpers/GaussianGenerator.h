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
  template <class Particle, class ParticleCell>
  static void fillWithParticles(autopas::AutoPas<Particle, ParticleCell> &autoPas, size_t numParticles,
                                Particle &defaultParicle = autopas::Particle(), double distributionMean = 5.0,
                                double distributionStdDev = 2.0) {
    std::default_random_engine generator(42);
    std::normal_distribution<double> distribution(distributionMean, distributionStdDev);

    for (size_t id = 0; id < numParticles;) {
      std::array<double, 3> position = {distribution(generator), distribution(generator), distribution(generator)};
      // only increment loop var (and place particle) if position is valid
      if (not autopas::utils::inBox(position, {0, 0, 0}, autoPas.getBoxMax())) continue;
      Particle p(defaultParicle);
      p.setR(position);
      p.setID(id++);
      autoPas.addParticle(p);
    }
  }
};
