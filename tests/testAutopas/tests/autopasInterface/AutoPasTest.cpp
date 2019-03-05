/**
 * @file AutoPasTest.cpp
 * @author seckler
 * @date 29.05.18
 */

#include "AutoPasTest.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Return;

/**
 * test whether the RegionIterator can be generated and something is returned.
 * This mainly makes certain, that the specific parts of the code can be compiled.
 */
TEST_F(AutoPasTest, getRegionParticleIterator) {
  auto iter = autoPas.getRegionIterator({0., 0., 0.}, {4., 4., 4.});

  for (; iter.isValid(); ++iter) {
    iter->setR(iter->getR());
  }
}

/**
 * Check whether an AutoPas object can be rebuild using the following strategy:
 * 0. create normal AutoPas container + initialize and fill with particles
 * 1. create additional AutoPas container + initialize
 * 2. fill additional AutoPas container with particles
 * 3. std::move the new container over the old one.
 */
TEST_F(AutoPasTest, checkRebuildingNewMove) {
  auto boxMax1 = autoPas.getBoxMax();
  ASSERT_DOUBLE_EQ(boxMax1[0], 10.);
  ASSERT_DOUBLE_EQ(boxMax1[1], 10.);
  ASSERT_DOUBLE_EQ(boxMax1[2], 10.);
  // 0. add some particles
  Particle p1({1., 1., 1.}, {0., 0., 0.}, 0);
  autoPas.addParticle(p1);
  Particle p2({2., 2., 2.}, {0., 0., 0.}, 1);
  autoPas.addParticle(p2);
  {
    // 1. create new AutoPas container + initialize
    decltype(autoPas) autoPasTmp;
    autoPasTmp.init({0., 0., 0.}, {5., 5., 5.}, 1., 0, 1, {autopas::ContainerOptions::linkedCells},
                    {autopas::TraversalOptions::c08});

    // ensure no particles
    for (auto iter = autoPasTmp.begin(); iter.isValid(); ++iter) {
      FAIL() << "there shouldn't already be any particles in the container, but there are." << std::endl;
    }

    // 2. copy particles
    for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
      autoPasTmp.addParticle(*iter);
    }

    {
      auto iter = autoPasTmp.begin();
      ASSERT_TRUE(iter.isValid());
      auto id1 = iter->getID();
      ++iter;
      ASSERT_TRUE(iter.isValid());
      ASSERT_EQ(id1 + iter->getID(), 1);
    }

    // 3. move new container over old one. will loose information of the old one.
    autoPas = std::move(autoPasTmp);
  }

  // ensure same particles as before
  {
    auto iter = autoPas.begin();
    ASSERT_TRUE(iter.isValid());
    auto id1 = iter->getID();
    ++iter;
    ASSERT_TRUE(iter.isValid());
    ASSERT_EQ(id1 + iter->getID(), 1);
  }

  // ensure information has changed
  auto boxMax2 = autoPas.getBoxMax();
  ASSERT_DOUBLE_EQ(boxMax2[0], 5.);
  ASSERT_DOUBLE_EQ(boxMax2[1], 5.);
  ASSERT_DOUBLE_EQ(boxMax2[2], 5.);

  // ensure logger still working
  AutoPasLog(info, "test logger working.");
}

/**
 * Check whether an AutoPas object can be rebuild using the following strategy:
 * 0. create normal AutoPas container + initialize and fill with particles
 * 1. copy particles out of first container
 * 2. recreate container by calling constructor + initialize again
 * 3. fill container with copied particles
 */
TEST_F(AutoPasTest, checkRebuildingCopyCreateNew) {
  auto boxMax1 = autoPas.getBoxMax();
  ASSERT_DOUBLE_EQ(boxMax1[0], 10.);
  ASSERT_DOUBLE_EQ(boxMax1[1], 10.);
  ASSERT_DOUBLE_EQ(boxMax1[2], 10.);

  // 0. add some particles
  Particle p1({1., 1., 1.}, {0., 0., 0.}, 0);
  autoPas.addParticle(p1);
  Particle p2({2., 2., 2.}, {0., 0., 0.}, 1);
  autoPas.addParticle(p2);
  {
    std::vector<Particle> particleVector;
    particleVector.reserve(autoPas.getNumberOfParticles());

    // 1. copy particles out of first container
    for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
      particleVector.push_back(*iter);
    }
    // 2. recreate container by calling constructor + init
    autoPas = decltype(autoPas)();
    autoPas.init({0., 0., 0.}, {5., 5., 5.}, 1., 0, 1, {autopas::ContainerOptions::linkedCells},
                 {autopas::TraversalOptions::c08});

    // ensure no particles
    for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
      FAIL() << "there shouldn't already be any particles in the container, but there are." << std::endl;
    }

    // 3. copy particles from particleVector into container
    for (auto particle : particleVector) {
      autoPas.addParticle(particle);
    }
  }

  // ensure same particles as before
  {
    auto iter = autoPas.begin();
    ASSERT_TRUE(iter.isValid());
    auto id1 = iter->getID();
    ++iter;
    ASSERT_TRUE(iter.isValid());
    ASSERT_EQ(id1 + iter->getID(), 1);
  }

  // ensure information has changed
  auto boxMax2 = autoPas.getBoxMax();
  ASSERT_DOUBLE_EQ(boxMax2[0], 5.);
  ASSERT_DOUBLE_EQ(boxMax2[1], 5.);
  ASSERT_DOUBLE_EQ(boxMax2[2], 5.);

  // ensure logger still working
  AutoPasLog(info, "test logger working.");
}

TEST_F(AutoPasTest, checkNeedsContainerUpdate) {
  // for linked cells this should be false
  EXPECT_TRUE(autoPas.needsContainerUpdate());

  // now build verlet lists
  autoPas.init({0., 0., 0.}, {5., 5., 5.}, 1., 0, 2, {autopas::ContainerOptions::verletLists}, {});
  // after build this should be false
  EXPECT_TRUE(autoPas.needsContainerUpdate());

  // run once, builds verlet lists. (here for 0 particles)
  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, allowsNewton3()).Times(AtLeast(1)).WillOnce(Return(true));
  EXPECT_CALL(emptyFunctor, allowsNonNewton3()).Times(AtLeast(1)).WillOnce(Return(false));
  EXPECT_CALL(emptyFunctor, isRelevantForTuning()).Times(AtLeast(1)).WillOnce(Return(true));
  autopas::C08Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> dummyTraversal({0, 0, 0},
                                                                                               &emptyFunctor);
  autoPas.iteratePairwise(&emptyFunctor, autopas::aos);

  // now verlet lists should be valid.
  EXPECT_FALSE(autoPas.needsContainerUpdate());
}
