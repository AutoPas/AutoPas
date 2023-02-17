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
  Molecule p1({1., 1., 1.}, {0., 0., 0.}, 0);
  autoPas.addParticle(p1);
  Molecule p2({2., 2., 2.}, {0., 0., 0.}, 1);
  autoPas.addParticle(p2);
  {
    // 1. create new AutoPas container + initialize
    decltype(autoPas) autoPasTmp;
    autoPasTmp.setBoxMin({0., 0., 0.});
    autoPasTmp.setBoxMax({5., 5., 5.});
    autoPasTmp.setCutoff(1.);
    autoPasTmp.setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPasTmp.setAllowedTraversals({autopas::TraversalOption::lc_c08});
    autoPasTmp.setOutputSuffix("tmp_");
    autoPasTmp.init();

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
  AutoPasLog(INFO, "test logger working.");
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
  Molecule p1({1., 1., 1.}, {0., 0., 0.}, 0);
  autoPas.addParticle(p1);
  Molecule p2({2., 2., 2.}, {0., 0., 0.}, 1);
  autoPas.addParticle(p2);
  {
    std::vector<Molecule> particleVector;
    particleVector.reserve(autoPas.getNumberOfParticles());

    // 1. copy particles out of first container
    for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
      particleVector.push_back(*iter);
    }
    // 2. recreate container by calling constructor + init
    autoPas = decltype(autoPas)();

    autoPas.setBoxMin({0., 0., 0.});
    autoPas.setBoxMax({5., 5., 5.});
    autoPas.setCutoff(1.);
    autoPas.setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPas.setAllowedTraversals({autopas::TraversalOption::lc_c08});
    autoPas.init();

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
  AutoPasLog(INFO, "test logger working.");
}

TEST_F(AutoPasTest, checkArgumentValidation) {
  // cutoff
  EXPECT_ANY_THROW(autoPas.setCutoff(-5.0));
  EXPECT_ANY_THROW(autoPas.setCutoff(0.0));
  EXPECT_NO_THROW(autoPas.setCutoff(0.5));

  // cell size
  EXPECT_ANY_THROW(autoPas.setCellSizeFactor(-5.0));
  EXPECT_ANY_THROW(autoPas.setCellSizeFactor(0.0));
  EXPECT_NO_THROW(autoPas.setCellSizeFactor(0.5));
}

template <typename AP>
void testConstIterator(const AP &ap, int numParticles) {
  int num = 0;
  for (auto iter = ap.begin(); iter != ap.end(); ++iter) {
    ++num;
  }
  EXPECT_EQ(num, numParticles);
}

TEST_F(AutoPasTest, checkConstIterator) {
  // with 0 particles
  testConstIterator(autoPas, autoPas.getNumberOfParticles());

  Molecule p1({1., 1., 1.}, {0., 0., 0.}, 0);
  autoPas.addParticle(p1);
  Molecule p2({2., 2., 2.}, {0., 0., 0.}, 1);
  autoPas.addParticle(p2);
  // with 2 particles
  testConstIterator(autoPas, autoPas.getNumberOfParticles());
}

void AutoPasTest::expectedParticles(size_t expectedOwned, size_t expectedHalo) {
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), expectedHalo + expectedOwned);
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned), expectedOwned);
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::halo), expectedHalo);
}

TEST_F(AutoPasTest, getNumParticlesTest) {
  // there should be no particles in an empty container
  expectedParticles(0, 0);

  Molecule particle{};

  // add a particle in the domain -> owned
  particle.setR({1, 1, 1});
  autoPas.addParticle(particle);
  expectedParticles(1, 0);

  // add a particle outside the domain -> halo
  particle.setR({-0.1, -0.1, -0.1});
  autoPas.addHaloParticle(particle);
  expectedParticles(1, 1);

  // update container is expected to delete all halo particles and return leaving particles
  auto leavingParticles = autoPas.updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0);
  expectedParticles(1, 0);

  // move the owned particle in the halo
  autoPas.begin()->setR({-0.2, -0.2, -0.2});
  leavingParticles = autoPas.updateContainer();
  EXPECT_EQ(leavingParticles.size(), 1);
  expectedParticles(0, 0);
}

TEST_F(AutoPasTest, getNumParticlesIteratorTest) {
  // there should be no particles in an empty container
  expectedParticles(0, 0);

  Molecule particle{};

  // add a particle in the domain -> owned
  int numParticles = 0;
  for (; numParticles < 5; ++numParticles) {
    particle.setR({(double)numParticles, (double)numParticles, (double)numParticles});
    autoPas.addParticle(particle);
    expectedParticles(numParticles + 1, 0);
  }

  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    autoPas.deleteParticle(iter);
    --numParticles;
    expectedParticles(numParticles, 0);
  }
}