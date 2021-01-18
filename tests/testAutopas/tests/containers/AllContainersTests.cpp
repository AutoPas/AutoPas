/**
 * @file TestsAllContainers.cpp
 * @author humig
 * @date 08.07.2019
 */

#include "AllContainersTests.h"

#include <gtest/gtest.h>

INSTANTIATE_TEST_SUITE_P(Generated, AllContainersTests, testing::ValuesIn(autopas::ContainerOption::getAllOptions()),
                         AllContainersTests::getParamToStringFunction());

/**
 * Checks if ParticleContainerInterface::getNumParticle() returns the correct number of particles.
 */
TEST_P(AllContainersTests, testGetNumParticles) {
  auto container = getInitializedContainer();
  EXPECT_EQ(container->getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  container->addParticle(p);
  EXPECT_EQ(container->getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  container->addParticle(p2);
  EXPECT_EQ(container->getNumParticles(), 2);
}

/**
 * Checks if ParticleContainerInterface::deleteAllParticles() deletes all particles.
 */
TEST_P(AllContainersTests, testDeleteAllParticles) {
  auto container = getInitializedContainer();
  EXPECT_EQ(container->getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  container->addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  container->addParticle(p2);
  EXPECT_EQ(container->getNumParticles(), 2);

  container->deleteAllParticles();
  EXPECT_EQ(container->getNumParticles(), 0);
}

/**
 * Checks if addParticle() only accepts particles in the domain and throws otherwise, and if addHaloParticle() never
 * throws.
 */
TEST_P(AllContainersTests, testParticleAdding) {
  auto container = getInitializedContainer();
  int id = 1;
  for (double x : {boxMin[0] - 1.5, boxMin[0] - .5, boxMin[0], boxMin[0] + 5., boxMax[0] - 0.001, boxMax[0],
                   boxMax[0] + .5, boxMax[0] + 1.5}) {
    for (double y : {boxMin[1] - 1.5, boxMin[1] - .5, boxMin[1], boxMin[1] + 5., boxMax[1] - 0.001, boxMax[1],
                     boxMax[1] + .5, boxMax[1] + 1.5}) {
      for (double z : {boxMin[2] - 1.5, boxMin[2] - .5, boxMin[2], boxMin[2] + 5., boxMax[2] - 0.001, boxMax[2],
                       boxMax[2] + .5, boxMax[2] + 1.5}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        if (x == -1.5 or y == -1.5 or z == -1.5 or x == 11.5 or y == 11.5 or z == 11.5) {
          EXPECT_ANY_THROW(container->addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(container->addHaloParticle(p));  // much outside, still ok because it is ignored!
        } else if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
                   z == 10.5) {
          EXPECT_ANY_THROW(container->addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(container->addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(container->addParticle(p));  // inside, therefore ok!
          EXPECT_ANY_THROW(
              // inside, and not ok, as halo particles cannot be added inside of the domain!
              container->addHaloParticle(p));
        }
      }
    }
  }
}

/**
 * Checks if updateContainer() deletes particles in halo.
 */
TEST_P(AllContainersTests, testUpdateContainerHalo) {
  auto container = getInitializedContainer();
  autopas::Particle p({boxMin[0] - 0.5, boxMin[1] - 0.5, boxMin[2] - 0.5}, {0, 0, 0}, 42);
  container->addHaloParticle(p);

  EXPECT_EQ(container->getNumParticles(), 1);
  EXPECT_EQ(container->begin()->getID(), 42);

  auto invalidParticles = container->updateContainer();

  // no particle should be returned
  EXPECT_EQ(invalidParticles.size(), 0);

  // no particle should remain
  auto iter = container->begin();
  EXPECT_FALSE(iter.isValid());
}

/**
 * Checks if updateContainer deletes dummy particles.
 * @param previouslyOwned Specifies whether the particle was previously owned.
 */
void AllContainersTests::testUpdateContainerDeletesDummy(bool previouslyOwned) {
  static std::atomic<unsigned long> numParticles = 0;

  class TestParticle : public autopas::Particle {
   public:
    TestParticle(std::array<double, 3> r, std::array<double, 3> v, unsigned long id) : Particle(r, v, id) {
      ++numParticles;
    }
    TestParticle(const TestParticle &testParticle) : Particle(testParticle) { ++numParticles; }
    ~TestParticle() override { --numParticles; }
  };

  // We need the container to use TestParticle!
  auto container = getInitializedContainer<TestParticle>();

  // Add particle
  {
    double pos = previouslyOwned ? 0.5 : -0.5;
    TestParticle p({pos, pos, pos}, {0, 0, 0}, 42);
    if (previouslyOwned) {
      container->addParticle(p);
    } else {
      container->addHaloParticle(p);
    }
  }
  // Mark particle as deleted
  {
    auto iter = container->begin();
    ASSERT_TRUE(iter.isValid());
    autopas::internal::markParticleAsDeleted(*iter);
  }
  // Check that we do not iterate over it.
  {
    auto iter = container->begin();
    ASSERT_FALSE(iter.isValid());
  }

  // This should remove the dummy particle(s), while not returning it as invalid particle.
  auto invalidParticles = container->updateContainer();

  // No particle should be returned
  EXPECT_EQ(invalidParticles.size(), 0);
  // The particle should no longer exist, as it should be cleared.
  EXPECT_EQ(numParticles, 0);

  // no particle should remain, therefore the iterator should be invalid!
  EXPECT_FALSE(container->begin().isValid());

  container->deleteAllParticles();
  ASSERT_EQ(numParticles, 0);
}

/**
 * Checks if updateContainer deletes dummy particles.
 */
TEST_P(AllContainersTests, testUpdateContainerDeletesPreviouslyOwnedDummy) { testUpdateContainerDeletesDummy(true); }

TEST_P(AllContainersTests, testUpdateContainerDeletesPreviouslyHaloDummy) { testUpdateContainerDeletesDummy(false); }