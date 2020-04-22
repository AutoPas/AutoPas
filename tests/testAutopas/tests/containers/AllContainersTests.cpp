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
 * Checks if ParticleContainer::getNumParticle() returns the correct number of particles.
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
          EXPECT_ANY_THROW(
              // inside, and not ok, as halo particles cannot be added inside of the domain!
              _container->addHaloParticle(p));
        }
      }
    }
  }
}

/**
 * Checks if updateContainer() deletes particles in halo.
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