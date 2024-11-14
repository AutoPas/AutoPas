
#include "tests/zonalMethods/HalfShellTest.h"

#include <gtest/gtest.h>

void HalfShellTest::initContainer(AutoPasType &autopas, std::vector<ParticleType> particles) {
  autopas.setBoxMin(_boxMin);
  autopas.setBoxMax(_boxMax);
  autopas.setCutoff(_cutoff);
  autopas.init();

  // insert particles
  for (auto &particle : particles) {
    autopas.addParticle(particle);
  }
}

/**
 * Check if the import and export regions got initialized accordingly,
 * also if the collecting the particles yields the expected results
 */
TEST_F(HalfShellTest, regionInitialization) {
  ASSERT_EQ(_exportRegions.size(), HalfShell::_regionCount);
  ASSERT_EQ(_importRegions.size(), HalfShell::_regionCount);

  // homogenously distribute particles with a distance of 1 into the 10x10x10 domain
  std::vector<ParticleType> particles;
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 10; ++j) {
      for (int k = 0; k < 10; ++k) {
        ParticleType particle;
        particle.setR({static_cast<double>(i) + 0.5, static_cast<double>(j) + 0.5, static_cast<double>(k) + 0.5});
        particles.push_back(particle);
      }
    }
  }

  // add particles to container
  initContainer(_autopas, particles);

  // collect particles for each export region
  collectParticles(_autopas);

  // NOTE: this is a relatively weak check, but it suffices for now.
  size_t particleCount = 0;
  for (auto buffer : _regionBuffers) {
    particleCount += buffer.size();
  }

  /*
   *  Calculation:
   *  # z > 0 regions
   *  Region( origin: [0, 0, 9] | size: [10, 10, 1] ) -> 100
   *  Region( origin: [0, 0, 9] | size: [10, 1, 1] ) -> 10
   *  Region( origin: [0, 0, 9] | size: [1, 10, 1] ) -> 10
   *  Region( origin: [0, 0, 9] | size: [1, 1, 1] ) -> 1
   *  Region( origin: [0, 9, 9] | size: [10, 1, 1] ) -> 10
   *  Region( origin: [0, 9, 9] | size: [1, 1, 1] ) -> 1
   *  Region( origin: [9, 0, 9] | size: [1, 10, 1] ) -> 10
   *  Region( origin: [9, 0, 9] | size: [1, 1, 1] ) -> 1
   *  Region( origin: [9, 9, 9] | size: [1, 1, 1] ) -> 1
   *  # z = 0 and y > 0 regions
   *  Region( origin: [0, 9, 0] | size: [10, 1, 10] ) -> 100
   *  Region( origin: [0, 9, 0] | size: [1, 1, 10] ) -> 10
   *  Region( origin: [9, 9, 0] | size: [1, 1, 10] ) -> 10
   *  #  z = 0 and y = 0 and x > 0 regions
   *  Region( origin: [9, 0, 0] | size: [1, 10, 10] ) -> 100
   *  Sum of particles:
   *  100 + 10 + 10 + 1 + 100 + 10 + 10 + 1 + 100 + 10 + 1 + 10 + 1 = 364
   */
  ASSERT_EQ(particleCount, 364);
}
