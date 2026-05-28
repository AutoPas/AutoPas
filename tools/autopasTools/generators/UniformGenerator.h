/**
 * @file UniformGenerator.h
 * @author seckler
 * @date 22.05.18
 */

#pragma once

#include <array>
#include <random>

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/inBox.h"

namespace autopasTools::generators {
/**
 * Generator class for uniform distributions
 */
namespace UniformGenerator {
/**
 * Generate a random position within a given box.
 * @param generator Random engine
 * @param boxMin
 * @param boxMax
 * @return 3D array with random values
 */
std::array<double, 3> randomPosition(std::mt19937 &generator, const std::array<double, 3> &boxMin,
                                     const std::array<double, 3> &boxMax);

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
void fillWithParticles(Container &container, const Particle &defaultParticle, const std::array<double, 3> &boxMin,
                       const std::array<double, 3> &boxMax, unsigned long numParticles = 100ul, unsigned int seed = 42);

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
void fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                           unsigned long numParticles, const HaloAddFunction &haloAddFunction, unsigned int seed = 42);

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
void fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                           unsigned long numParticles, unsigned int seed = 42) {
  fillWithHaloParticles(
      container, defaultParticle, haloWidth, numParticles,
      [](Container &c, Particle &p) {
        p.setOwnershipState(autopas::OwnershipState::halo);
        c.addHaloParticle(p);
      },
      seed);
}
};  // namespace UniformGenerator

template <class Container, class Particle>
void UniformGenerator::fillWithParticles(Container &container, const Particle &defaultParticle,
                                         const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                                         unsigned long numParticles, unsigned int seed) {
  std::mt19937 generator(seed);

  for (unsigned long i = defaultParticle.getID(); i < defaultParticle.getID() + numParticles; ++i) {
    Particle particle(defaultParticle);
    auto pos = randomPosition(generator, boxMin, boxMax);
    particle.setR({
      static_cast<Particle::ParticleSoAFloatPrecision>(pos.at(0)),
      static_cast<Particle::ParticleSoAFloatPrecision>(pos.at(1)),
      static_cast<Particle::ParticleSoAFloatPrecision>(pos.at(2))
    });
    particle.setID(i);
    particle.setOwnershipState(autopas::OwnershipState::owned);
    container.addParticle(particle);
  }
}

template <class Container, class Particle, class HaloAddFunction>
void UniformGenerator::fillWithHaloParticles(Container &container, const Particle &defaultParticle, double haloWidth,
                                             unsigned long numParticles, const HaloAddFunction &haloAddFunction,
                                             unsigned int seed) {
  std::mt19937 generator(seed);

  auto haloBoxMin = container.getBoxMin();
  auto haloBoxMax = container.getBoxMax();

  // increase the box size not exactly by the width to make it exclusive
  for (unsigned int i = 0; i < 3; ++i) {
    haloBoxMin[i] -= haloWidth * .99;
    haloBoxMax[i] += haloWidth * .99;
  }

  for (unsigned long i = defaultParticle.getID(); i < defaultParticle.getID() + numParticles; ++i) {
    auto pos = randomPosition(generator, haloBoxMin, haloBoxMax);
    // we only want to add particles not in the actual box
    while (autopas::utils::inBox(pos, container.getBoxMin(), container.getBoxMax())) {
      pos = randomPosition(generator, haloBoxMin, haloBoxMax);
    }
    Particle particle(defaultParticle);
    particle.setR(pos);
    particle.setID(i);
    haloAddFunction(container, particle);
  }
}
}  // namespace autopasTools::generators
