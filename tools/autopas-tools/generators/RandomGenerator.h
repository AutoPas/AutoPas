/**
 * @file RandomGenerator.h
 * @author seckler
 * @date 22.05.18
 */

#pragma once

#include <array>

#include "autopas/utils/inBox.h"

namespace autopas_tools::generators {
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
  static std::array<double, 3> randomPosition(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax);

  /**
   * Fills any container (also AutoPas object) with randomly uniformly distributed particles.
   * Particle properties will be used from the default particle. Particle IDs start from the default particle.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param defaultParticle
   * @param boxMin
   * @param boxMax
   * @param numParticles
   * @param seed
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container &container, const Particle &defaultParticle,
                                const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                unsigned long numParticles = 100ul, unsigned int seed = 42);

  /**
   * Fills only a given part of a container (also AutoPas object) with randomly uniformly distributed particles.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param defaultParticle
   * @param haloWidth
   * @param numParticles
   * @param seed
   */
  template <class Container, class Particle>
  static void fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                                    unsigned long numParticles = 100ul, unsigned int seed = 42);
};

template <class Container, class Particle>
void RandomGenerator::fillWithParticles(Container &container, const Particle &defaultParticle,
                                        const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                        unsigned long numParticles, unsigned int seed) {
  srand(seed);  // fixed seedpoint

  for (unsigned long i = defaultParticle.getID(); i < defaultParticle.getID() + numParticles; ++i) {
    Particle particle(defaultParticle);
    particle.setR(randomPosition(boxMin, boxMax));
    particle.setID(i);
    container.addParticle(particle);
  }
}

template <class Container, class Particle>
void RandomGenerator::fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                                            unsigned long numParticles, unsigned int seed) {
  srand(seed);  // fixed seedpoint

  auto haloBoxMin = container.getBoxMin();
  auto haloBoxMax = container.getBoxMax();

  // increase the box size not exactly by the width to make it exclusive
  for (unsigned int i = 0; i < 3; ++i) {
    haloBoxMin[i] -= haloWidth * .99;
    haloBoxMax[i] += haloWidth * .99;
  }

  for (unsigned long i = defaultParticle.getID(); i < defaultParticle.getID() + numParticles; ++i) {
    auto pos = randomPosition(haloBoxMax, haloBoxMax);
    // we only want to add particles not in the actual box
    while (autopas::utils::inBox(pos, container.getBoxMin(), container.getBoxMax())) {
      pos = randomPosition(haloBoxMax, haloBoxMax);
    }
    Particle particle(defaultParticle);
    particle.setR(pos);
    particle.setID(i);
    container.addHaloParticle(particle);
  }
}
}  // namespace autopas_tools::generators
