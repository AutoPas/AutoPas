/**
 * @file RandomGenerator.h
 * @author seckler
 * @date 22.05.18
 */


#pragma once

#include <array>


class RandomGenerator {
 private:
  static double fRand(double fMin, double fMax);

  static std::array<double, 3> randomPosition(const std::array<double, 3>& boxMin,
                                       const std::array<double, 3>& boxMax);
 public:
  template <class Container, class Particle>
  static void fillWithParticles(Container& container, Particle defaultParticle) {
    srand(42);  // fixed seedpoint
    int numParticles = 100;

    for (int i = 0; i < numParticles; ++i) {
      auto id = static_cast<unsigned long>(i);
      Particle particle = defaultParticle;
      particle.setR(randomPosition(container.getBoxMin(), container.getBoxMax()));
      particle.setID(id);
      container.addParticle(particle);
    }
  }
};



