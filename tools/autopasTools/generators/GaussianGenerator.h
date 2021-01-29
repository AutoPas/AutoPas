/**
 * @file GaussianGenerator.h
 * @author F. Gratl
 * @date 7/31/18
 */

#pragma once

#include <random>

#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ParticleTypeTrait.h"

namespace autopasTools::generators {
/**
 * Generator class for gaussian distributions
 */
class GaussianGenerator {
 public:
  /**
   * Fills any container (also AutoPas object) with randomly 3D gaussian distributed particles.
   * Particle properties will be used from the default particle. Particle IDs start from the default particle.
   * @tparam Container Arbitrary container class. Needs to support addParticle().
   * @param container
   * @param boxMin
   * @param boxMax
   * @param numParticles number of particles
   * @param defaultParticle inserted particle
   * @param distributionMean mean value / expected value. Mean is relative to boxMin.
   * @param distributionStdDev standard deviation
   */
  template <class Container>
  static void fillWithParticles(Container &container, const std::array<double, 3> &boxMin,
                                const std::array<double, 3> &boxMax, size_t numParticles,
                                const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle =
                                    typename autopas::utils::ParticleTypeTrait<Container>::value(),
                                const std::array<double, 3> &distributionMean = {5., 5., 5.},
                                const std::array<double, 3> &distributionStdDev = {2., 2., 2.});

 private:
  /**
   * Maximum number of attempts the random generator gets to find a valid position before considering the input to be
   * bad
   */
  constexpr static size_t _maxAttempts = 100;
};

template <class Container>
void GaussianGenerator::fillWithParticles(
    Container &container, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, size_t numParticles,
    const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle,
    const std::array<double, 3> &distributionMean, const std::array<double, 3> &distributionStdDev) {
  std::default_random_engine generator(42);
  std::array<std::normal_distribution<double>, 3> distributions = {
      std::normal_distribution<double>{distributionMean[0], distributionStdDev[0]},
      std::normal_distribution<double>{distributionMean[1], distributionStdDev[1]},
      std::normal_distribution<double>{distributionMean[2], distributionStdDev[2]}};

  for (uint64_t i = defaultParticle.getID(); i < defaultParticle.getID() + numParticles; ++i) {
    std::array<double, 3> position = {distributions[0](generator), distributions[1](generator),
                                      distributions[2](generator)};
    // verifies that position is valid
    for (size_t attempts = 1; attempts <= _maxAttempts and (not autopas::utils::inBox(position, boxMin, boxMax));
         ++attempts) {
      if (attempts == _maxAttempts) {
        std::ostringstream errormessage;
        errormessage << "GaussianGenerator::fillWithParticles(): Could not find a valid position for particle " << i
                     << "after" << _maxAttempts << "attempts. Check if your parameters make sense:" << std::endl
                     << "BoxMin       = " << autopas::utils::ArrayUtils::to_string(boxMin) << std::endl
                     << "BoxMax       = " << autopas::utils::ArrayUtils::to_string(boxMax) << std::endl
                     << "Gauss mean   = " << autopas::utils::ArrayUtils::to_string(distributionMean) << std::endl
                     << "Gauss stdDev = " << autopas::utils::ArrayUtils::to_string(distributionStdDev) << std::endl;
        throw std::runtime_error(errormessage.str());
      }
      position = {distributions[0](generator), distributions[1](generator), distributions[2](generator)};
    };
    auto p = defaultParticle;
    p.setR(position);
    p.setID(i);
    container.addParticle(p);
  }
}
}  // namespace autopasTools::generators
