/**
 * @file RandomGenerator.h
 * @author seckler
 * @date 22.05.18
 */

#pragma once

#include <array>
#include "autopas/utils/inBox.h"

class RandomGenerator {
 private:
  /**
   * returns a random number between fMin and fMax
   * @param fMin min. number (inclusive)
   * @param fMax max. number (inclusive)
   * @return double
   */
  static double fRand(double fMin, double fMax);

 public:
  /**
   * returns a random position in the space between boxMin and boxMax
   * @param boxMin min. position
   * @param boxMax max. position
   * @return std::array<double, 3>
   */
  static std::array<double, 3> randomPosition(const std::array<double, 3>& boxMin, const std::array<double, 3>& boxMax);

  /**
   * fills the given container with random distributed particles
   * @tparam Container
   * @tparam Particle Type of particle to be generated
   * @param container
   * @param defaultParticle inserted particle
   * @param numParticles number of particles
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container& container, const Particle& defaultParticle,
                                unsigned long numParticles = 100ul);

  /**
   * fills the given container with random distributed particles between boxMin and boxMax
   * @tparam Container
   * @tparam Particle Type of particle to be generated
   * @param container
   * @param defaultParticle inserted particle
   * @param boxMin min. position
   * @param boxMax max. position
   * @param numParticles number of particles
   */
  template <class Container, class Particle>
  static void fillWithParticles(Container& container, const Particle& defaultParticle,
                                const std::array<double, 3>& boxMin, const std::array<double, 3>& boxMax,
                                unsigned long numParticles = 100ul);

  /**
   * fills halo of given container with n = numParticles particles
   * @tparam Container
   * @tparam Particle Type of particle to be generated
   * @param container
   * @param defaultParticle inserted particle
   * @param haloWidth width of halo
   * @param numParticles number of particles
   */
  template <class Container, class Particle>
  static void fillWithHaloParticles(Container& container, const Particle& defaultParticle, double haloWidth,
                                    unsigned long numParticles = 100ul);
};

template <class Container, class Particle>
void RandomGenerator::fillWithParticles(Container& container, const Particle& defaultParticle,
                                        unsigned long numParticles) {
  RandomGenerator::fillWithParticles(container, defaultParticle, container.getBoxMin(), container.getBoxMax(),
                                     numParticles);
}

template <class Container, class Particle>
void RandomGenerator::fillWithParticles(Container& container, const Particle& defaultParticle,
                                        const std::array<double, 3>& boxMin, const std::array<double, 3>& boxMax,
                                        unsigned long numParticles) {
  srand(42);  // fixed seedpoint

  for (unsigned long i = 0; i < numParticles; ++i) {
    Particle particle(defaultParticle);
    particle.setR(randomPosition(boxMin, boxMax));
    particle.setID(i);
    container.addParticle(particle);
  }
}

template <class Container, class Particle>
void RandomGenerator::fillWithHaloParticles(Container& container, const Particle& defaultParticle, double haloWidth,
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
