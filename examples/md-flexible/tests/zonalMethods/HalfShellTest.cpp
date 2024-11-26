
#include "tests/zonalMethods/HalfShellTest.h"

#include <gtest/gtest.h>

#include "src/TypeDefinitions.h"

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

/**
 * Test if the result is recollected correctly
 * NOTE: this test assumes that the region initialization test passes
 */
TEST_F(HalfShellTest, testResultRecollection) {
  initContainer(_autopas, {});

  auto hsCondition = [](const int d[3]) {
    /**
     * Stencil:
     *  z > 0 +
     *  z == 0 and y > 0 +
     *  z == 0 and y == 0 and x > 0
     */
    return d[2] > 0 or (d[2] == 0 and (d[1] > 0 or (d[1] == 0 and d[0] > 0)));
  };

  std::map<std::tuple<int, int, int>, std::vector<size_t>> regionParticles;
  size_t id = 0;
  // iterate over neighbours
  for (int x = -1; x < 2; ++x) {
    for (int y = -1; y < 2; ++y) {
      for (int z = -1; z < 2; ++z) {
        int d[3] = {x, y, z};
        // if neighbour is in the stencil
        if (hsCondition(d)) {
          std::vector<ParticleType> particles;
          std::vector<size_t> indices;
          // set a random number of particles
          size_t numberOfParticles = rand() % 100 + 1;
          for (size_t i = 0; i < numberOfParticles; i++) {
            double px = (x == -1) * (_boxMin[0] - 0.5) + (x == 0) * (_boxMin[0] + 5) + (x == 1) * (_boxMax[0] + 0.5);
            double py = (y == -1) * (_boxMin[0] - 0.5) + (y == 0) * (_boxMin[0] + 5) + (y == 1) * (_boxMax[0] + 0.5);
            double pz = (z == -1) * (_boxMin[0] - 0.5) + (z == 0) * (_boxMin[0] + 5) + (z == 1) * (_boxMax[0] + 0.5);
            ParticleType p({px, py, pz}, {0, 0, 0}, id);
            indices.push_back(id++);
            particles.push_back(p);
          }
          _autopas.addHaloParticles(particles);
          regionParticles.insert_or_assign(std::make_tuple(x, y, z), indices);
        }
      }
    }
  }

  // recollect
  recollectResultsFromContainer(_autopas);

  // check if particles are in their correct region
  size_t index = 0;
  for (auto &imRegion : _importRegions) {
    auto expectedParticles = regionParticles.at(
        std::make_tuple(imRegion.getNeighbour()[0], imRegion.getNeighbour()[1], imRegion.getNeighbour()[2]));
    for (auto &p : _importBuffers[index]) {
      bool found = false;
      size_t r_index = 0;
      for (auto &r : expectedParticles) {
        if (p.getID() == r) {
          found = true;
          expectedParticles.erase(expectedParticles.begin() + r_index);
          break;
        }
        ++r_index;
      }
      EXPECT_TRUE(found);
    }
    EXPECT_EQ(expectedParticles.size(), 0);
    ++index;
  }
}
