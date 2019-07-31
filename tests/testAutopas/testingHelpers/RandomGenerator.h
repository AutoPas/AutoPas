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
  static std::array<double, 3> randomPosition(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax);

  /**
   * Fills any container (also AutoPas object) with randomly uniformly distributed particles.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param defaultParticle
   * @param numParticles
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container &container, const Particle &defaultParticle,
                                unsigned long numParticles = 100ul,
                                const std::array<double, 3> &velocity = {0., 0., 0.});

  /**
   * Fills the given container with randomly distributed particles between boxMin and boxMax.
   * @tparam Container
   * @tparam Particle Type of particle to be generated
   * @param container
   * @param typeId
   * @param id
   * @param defaultParticle inserted particle
   * @param boxMin min. position
   * @param boxMax max. position
   * @param numParticles number of particles
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container &container,size_t typeId,size_t id,const Particle &defaultParticle,
                                const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                unsigned long numParticles = 100ul,
                                const std::array<double, 3> &velocity = {0., 0., 0.});

  /**
   * Fills only a given part of a container (also AutoPas object) with randomly uniformly distributed particles.
   * @tparam Container Arbitrary container class that needs to support getBoxMax() and addParticle().
   * @tparam Particle Type of the default particle.
   * @param container
   * @param defaultParticle
   * @param haloWidth
   * @param numParticles
   */
  template <class Container, class Particle>
  static void fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                                    unsigned long numParticles = 100ul);
};

template <class Container, class Particle>
void RandomGenerator::fillWithParticles(Container &container, const Particle &defaultParticle,
                                        unsigned long numParticles, const std::array<double, 3> &velocity) {
  RandomGenerator::fillWithParticles(container, defaultParticle, container.getBoxMin(), container.getBoxMax(),
                                     numParticles);
}

template <class Container, class Particle>
void RandomGenerator::fillWithParticles(Container &container,size_t typeId,size_t id, const Particle &defaultParticle,
                                        const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                        unsigned long numParticles, const std::array<double, 3> &velocity) {
  srand(42);  // fixed seedpoint

  for (unsigned long i = 0; i < numParticles; ++i) {
    Particle particle(defaultParticle);
    particle.setR(randomPosition(boxMin, boxMax));
    particle.setID(id);
    particle.setTypeId(typeId);
    particle.setV(velocity);
    container.addParticle(particle);
    id++;
  }
}

template <class Container, class Particle>
void RandomGenerator::fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                                            unsigned long numParticles) {
  srand(42);  // fixed seedpoint

  auto haloBoxMin = container.getBoxMin();
  auto haloBoxMax = container.getBoxMax();

  // increase the box size not exactly by the width to make it exclusive
  for (unsigned int i = 0; i < 3; ++i) {
    haloBoxMin[i] -= haloWidth * .99;
    haloBoxMax[i] += haloWidth * .99;
  }

  for (unsigned long i = 0; i < numParticles; ++i) {
    const auto pos = randomPosition(haloBoxMax, haloBoxMax);
    // we only want  to add particles not in the actual box
    if (autopas::utils::inBox(pos, container.getBoxMin(), container.getBoxMax())) continue;
    Particle particle(defaultParticle);
    particle.setR(pos);
    particle.setID(i);
    container.addHaloParticle(particle);
  }
}
