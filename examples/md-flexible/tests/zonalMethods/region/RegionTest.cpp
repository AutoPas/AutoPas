
#include "tests/zonalMethods/region/RegionTest.h"

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "src/zonalMethods/region/RectRegion.h"

void RegionTest::initContainer(AutoPasType &autopas, std::vector<ParticleType> particles) {
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({10., 10., 10.});
  autopas.setCutoff(1.);
  autopas.init();

  // insert particles
  for (auto &particle : particles) {
    autopas.addParticle(particle);
  }
}

/**
 * Test if RectRegion collects all particles in the container
 */
TEST_F(RegionTest, RectRegionCollectAllTest) {
  using namespace autopas::utils::ArrayMath;
  using namespace autopas::utils::ArrayUtils;

  // init particles to be inserted
  std::vector<ParticleType> particles;
  for (int i = 0; i < 5; ++i) {
    ParticleType particle;
    particle.setR({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
    particles.push_back(particle);
  }

  // insert particles and init container
  initContainer(_autopas, particles);

  // initialize RectRegion
  auto rectZone = RectRegion({0., 0., 0.}, {5., 5., 5.});

  // collect particles
  std::vector<ParticleType> buffer;
  rectZone.collectParticles(_autopas, buffer);

  unsigned long expectedSize = 5;
  // check that all particles were collected
  ASSERT_EQ(buffer.size(), expectedSize);
  for (int i = 0; i < expectedSize; ++i) {
    bool foundMatch = false;
    // find matching particle in particles
    for (int j = 0; j < particles.size(); ++j) {
      if (isNearRel(particles[j].getR(), buffer[i].getR())) {
        // remove from particles
        particles.erase(particles.begin() + j);
        foundMatch = true;
        break;
      }
    }
    ASSERT_TRUE(foundMatch) << "Particle collected by Zone not matching expected results: "
                            << to_string(buffer[i].getR()) << " not matched.";
  }
  // check if all particles were collected
  ASSERT_EQ(particles.size(), 0);
}

/**
 * Test if RectRegion collects some particles in the container
 * NOTE: this test showed that somehow the particle at {4, 4, 4} did not get collected,
 * despite being at the edge of the specified region.
 */
TEST_F(RegionTest, RectRegionCollectPartialTest) {
  using namespace autopas::utils::ArrayMath;
  using namespace autopas::utils::ArrayUtils;

  // init particles to be inserted
  std::vector<ParticleType> particles;
  for (int i = 0; i < 5; ++i) {
    ParticleType particle;
    particle.setR({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
    particles.push_back(particle);
  }

  // insert particles and init container
  initContainer(_autopas, particles);

  // initialize RectRegion
  auto rectZone = RectRegion({1., 0., 2.}, {3., 5., 2.});

  // collect particles
  std::vector<ParticleType> buffer;
  rectZone.collectParticles(_autopas, buffer);

  // fill in expected results
  particles.clear();
  for (int i = 2; i < 4; ++i) {
    ParticleType particle;
    particle.setR({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
    particles.push_back(particle);
  }

  unsigned long expectedSize = 2;

  // check that all particles were collected
  EXPECT_EQ(buffer.size(), expectedSize);
  for (int i = 0; i < expectedSize; ++i) {
    bool foundMatch = false;
    // find matching particle in particles
    for (int j = 0; j < particles.size(); ++j) {
      if (isNearRel(particles[j].getR(), buffer[i].getR())) {
        // remove from particles
        particles.erase(particles.begin() + j);
        foundMatch = true;
        break;
      }
    }
    ASSERT_TRUE(foundMatch) << "Particle collected by Zone not matching expected results: "
                            << to_string(buffer[i].getR()) << " not matched.";
  }
  // check if all particles were collected
  ASSERT_EQ(particles.size(), 0);
}
