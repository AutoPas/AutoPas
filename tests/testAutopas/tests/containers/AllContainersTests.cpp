/**
 * @file AllContainersTests.cpp
 * @author humig
 * @date 08.07.2019
 */

#include "AllContainersTests.h"

#include <gtest/gtest.h>

INSTANTIATE_TEST_SUITE_P(Generated, AllContainersTests,
                         ::testing::Combine(::testing::ValuesIn(autopas::ContainerOption::getAllOptions())),
                         AllContainersTests::getParamToStringFunction());

INSTANTIATE_TEST_SUITE_P(Generated, AllContainersTestsBothUpdates,
                         ::testing::Combine(::testing::ValuesIn(autopas::ContainerOption::getAllOptions()),
                                            ::testing::Bool()),
                         AllContainersTestsBothUpdates::getParamToStringFunction());

/**
 * Checks if ParticleContainerInterface::getNumParticle() returns the correct number of particles.
 */
TEST_P(AllContainersTests, testGetNumParticles) {
  auto container = getInitializedContainer(std::get<0>(GetParam()));
  EXPECT_EQ(container->getNumberOfParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  container->addParticle(p);
  EXPECT_EQ(container->getNumberOfParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  container->addParticle(p2);
  EXPECT_EQ(container->getNumberOfParticles(), 2);
}

/**
 * Checks if ParticleContainerInterface::deleteAllParticles() deletes all particles.
 */
TEST_P(AllContainersTests, testDeleteAllParticles) {
  auto container = this->getInitializedContainer(std::get<0>(GetParam()));
  EXPECT_EQ(container->getNumberOfParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  container->addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  container->addParticle(p2);
  EXPECT_EQ(container->getNumberOfParticles(), 2);

  container->deleteAllParticles();
  EXPECT_EQ(container->getNumberOfParticles(), 0);
}

/**
 * Checks if addParticle() only accepts particles in the domain and throws otherwise, and if addHaloParticle() never
 * throws.
 */
TEST_P(AllContainersTests, testParticleAdding) {
  auto container = getInitializedContainer(std::get<0>(GetParam()));
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
 * Add particles in the halo in every direction. Then delete them via deleteHaloParticles() and check that they are
 * gone.
 */
TEST_P(AllContainersTests, testDeleteHaloParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  auto container = getInitializedContainer(std::get<0>(GetParam()));

  std::array<double, 3> zeros{0, 0, 0};

  // counter and id for particles to be added
  size_t numParticles = 0;

  // calculate some distances needed later
  auto domainSize = container->getBoxMax() - container->getBoxMin();
  auto domainSizeHalf = domainSize * 0.5;
  auto domainCenter = container->getBoxMin() + domainSizeHalf;
  auto interactionLengthHalf = container->getInteractionLength() * 0.5;
  auto distCenterToMidHalo =
      domainSizeHalf + std::array<double, 3>{interactionLengthHalf, interactionLengthHalf, interactionLengthHalf};

  // setup: add some halo particles on every side and corner
  for (int x : {-1, 0, 1}) {
    for (int y : {-1, 0, 1}) {
      for (int z : {-1, 0, 1}) {
        if (x == 0 and y == 0 and z == 0) {
          continue;
        }
        auto pos = domainCenter + std::array<double, 3>{distCenterToMidHalo[0] * x, distCenterToMidHalo[1] * y,
                                                        distCenterToMidHalo[2] * z};
        Particle p{pos, zeros, numParticles++};
        container->addHaloParticle(p);
      }
    }
  }
  // sanity checks
  ASSERT_GT(numParticles, 0);
  ASSERT_EQ(container->getNumberOfParticles(), numParticles);

  // actual test:
  container->deleteHaloParticles();
  ASSERT_EQ(container->getNumberOfParticles(), 0);
}

/**
 * Checks if updateContainer() deletes particles in halo.
 */
TEST_P(AllContainersTestsBothUpdates, testUpdateContainerHalo) {
  auto container = getInitializedContainer(std::get<0>(GetParam()));
  autopas::Particle p({boxMin[0] - 0.5, boxMin[1] - 0.5, boxMin[2] - 0.5}, {0, 0, 0}, 42);
  container->addHaloParticle(p);

  EXPECT_EQ(container->getNumberOfParticles(), 1);
  EXPECT_EQ(container->begin()->getID(), 42);

  auto invalidParticles = container->updateContainer(std::get<1>(GetParam()));

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
void AllContainersTestsBothUpdates::testUpdateContainerDeletesDummy(bool previouslyOwned) {
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
  auto container = getInitializedContainer<TestParticle>(std::get<0>(GetParam()));

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
  auto invalidParticles = container->updateContainer(std::get<1>(GetParam()));

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
TEST_P(AllContainersTestsBothUpdates, testUpdateContainerDeletesPreviouslyOwnedDummy) {
  testUpdateContainerDeletesDummy(true);
}

TEST_P(AllContainersTestsBothUpdates, testUpdateContainerDeletesPreviouslyHaloDummy) {
  testUpdateContainerDeletesDummy(false);
}

/**
 * This test checks the correct behavior of the updateContainer(true) call.
 * Hereby we check, that:
 * a) Halo and leaving particles are converted to dummy particles and NOT completely removed.
 * b) Owned particles are not moved within the container, i.e., the pointer to the particle is still the same.
 */
TEST_P(AllContainersTests, testUpdateContainerKeepsNeighborListsValidIfSpecified) {
  auto container = getInitializedContainer(std::get<0>(GetParam()));

  {
    Particle p({-.1, -.1, -.1}, {0., 0., 0.}, 0);
    container->addHaloParticle(p);
  }

  {
    Particle p({.02, .1, .1}, {0., 0., 0.}, 1);
    container->addParticle(p);
  }

  {
    Particle p({1.23, .1, .1}, {0., 0., 0.}, 2);
    container->addParticle(p);
  }
  struct Values {
    unsigned long id;
    int occurrences;
  };
  std::map<Particle *, Values> previous_particles;
  for (auto &&p : *container) {
    // Iterates over owned and halo particles!
    previous_particles[&p] = {p.getID(), 0};
    if (p.getID() == 1) {
      // moves particle 1 outside the container (leaving particle!)
      p.addR({-0.04, 0., 0.});
    } else if (p.getID() == 2) {
      // Should normally move it to next linked cell, here the cell should NOT change!
      p.addR({0.04, 0., 0.});
    }
  }

  auto leavingParticles = container->updateContainer(true);

  // Particle 1 should be returned!
  ASSERT_EQ(leavingParticles.size(), 1ul);
  EXPECT_EQ(leavingParticles[0].getID(), 1);

  for (auto iter = container->begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
    ASSERT_EQ(previous_particles.count(&*iter), 1ul);

    ++previous_particles[&*iter].occurrences;
    if (iter->getID() == 0 or iter->getID() == 1) {
      // particle 0 and 1 should now be dummy!
      EXPECT_TRUE(iter->isDummy());
    } else {
      // particle 2 should still be owned!
      EXPECT_TRUE(iter->isOwned());
    }
  }

  // Every previous particle should be visited exactly once!
  for (auto [key, value] : previous_particles) {
    EXPECT_EQ(value.occurrences, 1) << "Particle with id " << value.id << " at address " << key
                                    << " not visited correctly!";
  }
}