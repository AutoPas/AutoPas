/**
 * @file RandomGenerator.h
 * @author seckler
 * @date 22.05.18
 */

#pragma once

#include <array>
#include "autopas/utils/inBox.h"

/**
 * Generator class for uniform distributions
 */
class RandomGenerator {
 private:
  /**
   * Simple random function
   * @param fMin
   * @param fMax
   * @return double between fMin and fMax
   */
  static double fRand(double fMin, double fMax);

 public:
  /**
   * Generate a random position within a given box.
   * @param boxMin
   * @param boxMax
   * @return 3D array with random values
   */
  static std::array<double, 3> randomPosition(const std::array<double, 3>& boxMin, const std::array<double, 3>& boxMax);

  /**
   * Fills any container (also AutoPas object) with randomly uniformly distributed particles.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param defaultParticle
   * @param numParticles
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container& container, Particle defaultParticle, unsigned long numParticles = 100) {
    srand(42);  // fixed seedpoint

    for (unsigned long i = 0; i < numParticles; ++i) {
      Particle particle = defaultParticle;
      particle.setR(randomPosition(container.getBoxMin(), container.getBoxMax()));
      particle.setID(i);
      container.addParticle(particle);
    }
  }

  /**
   * Fills only a given part of a container (also AutoPas object) with randomly uniformly distributed particles.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param defaultParticle
   * @param boxMin
   * @param boxMax
   * @param numParticles
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container& container, Particle defaultParticle, std::array<double, 3> boxMin,
                                std::array<double, 3> boxMax, int numParticles = 100) {
    srand(42);  // fixed seedpoint

    for (int i = 0; i < numParticles; ++i) {
      auto id = static_cast<unsigned long>(i);
      Particle particle = defaultParticle;
      particle.setR(randomPosition(boxMin, boxMax));
      particle.setID(id);
      container.addParticle(particle);
    }
  }

  /**
   * Fills the halo region of a container (also AutoPas object) with randomly uniformly distributed particles.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param defaultParticle
   * @param haloWidth
   * @param numParticles
   */
  template <class Container, class Particle>
  static void fillWithHaloParticles(Container& container, Particle defaultParticle, double haloWidth,
                                    unsigned long numParticles = 100) {
    srand(42);  // fixed seedpoint

    auto haloBoxMin = container.getBoxMin();
    auto haloBoxMax = container.getBoxMax();

    // increase the box size not exactly by the width to make it exclusive
    for (int i = 0; i < 3; ++i) {
      haloBoxMin[i] -= haloWidth * .99;
      haloBoxMax[i] += haloWidth * .99;
    }

    for (unsigned long i = 0; i < numParticles;) {
      Particle particle = defaultParticle;
      auto pos = randomPosition(haloBoxMax, haloBoxMax);
      // we only want  to add particles not in the actual box
      if (autopas::utils::inBox(pos, container.getBoxMin(), container.getBoxMax())) continue;
      particle.setR(pos);
      particle.setID(i);
      container.addHaloParticle(particle);
      ++i;
    }
  }
};
