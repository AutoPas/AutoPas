

#include "tests/zonalMethods/FullShellTest.h"

#include <gtest/gtest.h>

void FullShellTest::initContainer(AutoPasType &autopas, std::vector<ParticleType> particles) {
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
TEST_F(FullShellTest, testParticleCollection) {
  ASSERT_EQ(_exportRegions.size(), FullShell::_regionCount);
  ASSERT_EQ(_importRegions.size(), FullShell::_regionCount);

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
   *   #edge * 10 + #face * 100 + #corner * 1 =
   *   12 * 10 + 6 * 100 + 8 * 1 =
   *   120 + 600 + 8 = 728
   * */
  ASSERT_EQ(particleCount, 728);
}

