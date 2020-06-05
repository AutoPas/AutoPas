/**
 * @file RandomGenerator.h
 * @author seckler
 * @date 22.05.18
 */

#pragma once

#include <array>
#include <functional>
#include <iostream>

#include "autopas/utils/inBox.h"

namespace autopasTools::generators {
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

  /**
   * Detail class for fillWithHaloParticles
   * @tparam Container
   * @tparam Particle
   */
  template <class Container, class Particle>
  struct detail {
    /**
     * Default function to add halo particles to containers, uses function addHaloParticle().
     */
    constexpr static auto addHaloParticleF = [](Container &c, const Particle &p) { c.addHaloParticle(p); };
  };

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
   * Fills the halo of a container (also AutoPas object) with randomly uniformly distributed particles.
   * Use haloAddFunction to specify how to add particles to the container!
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @tparam HaloAddFunction function with signature void(Container&, Particle)
   * @param container
   * @param defaultParticle
   * @param haloWidth
   * @param numParticles
   * @param haloAddFunction the function of type HaloAddfunction that adds a halo particle to the container.
   * @param seed
   */
  template <class Container, class Particle, class HaloAddFunction>
  static void fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                                    unsigned long numParticles, const HaloAddFunction &haloAddFunction,
                                    unsigned int seed = 42);

  /**
   * Fills the halo of a container with randomly uniformly distributed particles.
   * Container needs to support addHaloParticle(). If it does not support this, please use the version above.
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
                                    unsigned long numParticles, unsigned int seed = 42) {
    fillWithHaloParticles(container, defaultParticle, haloWidth, numParticles,
                          detail<Container, Particle>::addHaloParticleF, seed);
  }
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

template <class Container, class Particle, class HaloAddFunction>
void RandomGenerator::fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                                            unsigned long numParticles, const HaloAddFunction &haloAddFunction,
                                            unsigned int seed) {
  srand(seed);  // fixed seedpoint

  auto haloBoxMin = container.getBoxMin();
  auto haloBoxMax = container.getBoxMax();

  // increase the box size not exactly by the width to make it exclusive
  for (unsigned int i = 0; i < 3; ++i) {
    haloBoxMin[i] -= haloWidth * .99;
    haloBoxMax[i] += haloWidth * .99;
  }

  for (unsigned long i = defaultParticle.getID(); i < defaultParticle.getID() + numParticles; ++i) {
    auto pos = randomPosition(haloBoxMin, haloBoxMax);
    // we only want to add particles not in the actual box
    while (autopas::utils::inBox(pos, container.getBoxMin(), container.getBoxMax())) {
      pos = randomPosition(haloBoxMin, haloBoxMax);
    }
    Particle particle(defaultParticle);
    particle.setR(pos);
    particle.setID(i);
    haloAddFunction(container, particle);
  }
}
}  // namespace autopasTools::generators
