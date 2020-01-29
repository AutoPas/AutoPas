/**
 * @file TestsAllContainers.cpp
 * @author humig
 * @date 08.07.2019
 */

#include "AllContainersTests.h"

#include <gtest/gtest.h>

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Invoke;
using ::testing::Values;

INSTANTIATE_TEST_SUITE_P(Generated, AllContainersTests, testing::ValuesIn([]() {
                           auto allOptions =
                               std::set<autopas::ContainerOption>{autopas::ContainerOption::getAllOptions()};
                           allOptions.erase(autopas::ContainerOption::verletClusterLists);
                           return allOptions;
                         }()),
                         AllContainersTests::getParamToStringFunction());

INSTANTIATE_TEST_SUITE_P(DISABLED_Generated, AllContainersTests,
                         testing::Values(autopas::ContainerOption::verletClusterLists),
                         AllContainersTests::getParamToStringFunction());

/**
 * Checks if ParticleContainer::getNumParticle() returns the correct number of particles.
 *
 * @Reviewer: Meant to replace
 *      VerletListsTest::testAddParticleNumParticle
 *      LinkedCellsTest::testGetNumParticles
 *      DirectSumContainerTest::testGetNumParticles
 */
TEST_P(AllContainersTests, testGetNumParticles) {
  EXPECT_EQ(_container->getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  _container->addParticle(p);
  EXPECT_EQ(_container->getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  _container->addParticle(p2);
  EXPECT_EQ(_container->getNumParticles(), 2);
}

/**
 * Checks if ParticleContainer::deleteAllParticles() deletes all particles.
 *
 * @Reviewer: Meant to replace
 * 	VerletListsTest::testDeleteAllParticles
 *      LinkedCellsTest::testDeleteAllParticles
 *      DirectSumContainerTest::testDeleteAllParticles
 */
TEST_P(AllContainersTests, testDeleteAllParticles) {
  EXPECT_EQ(_container->getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  _container->addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  _container->addParticle(p2);
  EXPECT_EQ(_container->getNumParticles(), 2);

  _container->deleteAllParticles();
  EXPECT_EQ(_container->getNumParticles(), 0);
}

/**
 * Checks if addParticle() only accepts particles in the domain and throws otherwise, and if addHaloParticle() never
 * throws.
 *
 * @Reviewer: Meant to replace
 *      LinkedCellsTest::testParticleAdding
 *      DirectSumContainerTest::testParticleAdding
 */
TEST_P(AllContainersTests, testParticleAdding) {
  int id = 1;
  for (double x : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
    for (double y : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
      for (double z : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        if (x == -1.5 or y == -1.5 or z == -1.5 or x == 11.5 or y == 11.5 or z == 11.5) {
          EXPECT_ANY_THROW(_container->addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(_container->addHaloParticle(p));  // much outside, still ok because it is ignored!
        } else if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
                   z == 10.5) {
          EXPECT_ANY_THROW(_container->addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(_container->addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(_container->addParticle(p));  // inside, therefore ok!
          EXPECT_NO_THROW(
              _container->addHaloParticle(p));  // inside, but still ok, as halo particle can be inside of the domain!
        }
      }
    }
  }
}

/**
 * Checks if the containers report that they need an update if a particle moved to halo.
 *
 * @Reviewer: Meant to partially replace (only the halo update, other conditions when it is necessary differ, right?)
 *      VerletListsTest::testIsContainerUpdateNeeded LinkedCellsTest::testIsContainerUpdateNeeded
 *	DirectSumContainerTest::testIsContainerUpdateNeeded
 */
TEST_P(AllContainersTests, testIsContainerUpdateNeededHalo) {
  EXPECT_FALSE(_container->isContainerUpdateNeeded());

  Particle p({1, 1, 1}, {0, 0, 0}, 0);
  _container->addParticle(p);
  EXPECT_FALSE(_container->isContainerUpdateNeeded());

  // Particle moves to halo cell -> needs update
  _container->begin()->setR({-1, -1, -1});
  EXPECT_TRUE(_container->isContainerUpdateNeeded());
}

/**
 * Checks if updateContainer() deletes particles in halo.
 *
 * @Reviewer: Meant to replace
 *      LinkedCellsTest::testUpdateContainerHalo
 *	DirectSumContainerTest::testUpdateContainerHalo
 */
TEST_P(AllContainersTests, testUpdateContainerHalo) {
  autopas::Particle p({-0.5, -0.5, -0.5}, {0, 0, 0}, 42);
  _container->addHaloParticle(p);

  EXPECT_EQ(_container->getNumParticles(), 1);
  EXPECT_EQ(_container->begin()->getID(), 42);

  auto invalidParticles = _container->updateContainer();

  // no particle should be returned
  EXPECT_EQ(invalidParticles.size(), 0);

  // no particle should remain
  auto iter = _container->begin();
  EXPECT_FALSE(iter.isValid());
}