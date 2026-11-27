/**
 * @file AutoPasTest.cpp
 * @author seckler
 * @date 29.05.18
 */

#include "AutoPasTest.h"

#include "autopas/AutoPasImpl.h"
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

TEST_F(AutoPasTest, getNumberOfParticlesTest) {
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

TEST_F(AutoPasTest, getNumberOfParticlesIteratorTest) {
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

/**
 * Verifies that bulk insertion methods (addParticles, addParticlesIf, addHaloParticles, addHaloParticlesIf)
 * accurately update the owned and halo particle counters in both rebuild mode (direct container insertion)
 * and non-rebuild mode (thread-private buffers and in-place halo updates).
 */
TEST_F(AutoPasTest, bulkAddParticlesCounterTest) {
  auto countByIterator = [&](const autopas::IteratorBehavior behavior) {
    size_t count = 0;
    for (auto iter = autoPas.begin(behavior); iter.isValid(); ++iter) {
      ++count;
    }
    return count;
  };

  auto checkCounts = [&](const size_t expectedOwned, const size_t expectedHalo) {
    expectedParticles(expectedOwned, expectedHalo);
    EXPECT_EQ(countByIterator(autopas::IteratorBehavior::owned), expectedOwned);
    EXPECT_EQ(countByIterator(autopas::IteratorBehavior::halo), expectedHalo);
  };

  // --- Phase 1: Rebuild state (iteration 0, neighbor lists invalid -> direct container insertion) ---
  // Advance to iteration 0 before the first computeInteractions so the next updateContainer enters iteration 1
  // (non-rebuild state with valid neighbor lists).
  EXPECT_TRUE(autoPas.updateContainer().empty());

  std::vector<Molecule> ownedBatch1;
  for (size_t i = 0; i < 20; ++i) {
    ownedBatch1.emplace_back(std::array<double, 3>{1.0 + 0.2 * static_cast<double>(i), 1.0, 1.0},
                             std::array<double, 3>{0.0, 0.0, 0.0}, i);
  }
  autoPas.addParticles(ownedBatch1);
  checkCounts(20, 0);

  std::vector<Molecule> ownedBatch2;
  for (size_t i = 20; i < 34; ++i) {
    ownedBatch2.emplace_back(std::array<double, 3>{1.0 + 0.2 * static_cast<double>(i - 20), 2.0, 1.0},
                             std::array<double, 3>{0.0, 0.0, 0.0}, i);
  }
  // Only even IDs (7 particles) pass predicate.
  autoPas.addParticlesIf(ownedBatch2, [](const auto &p) { return p.getID() % 2 == 0; });
  checkCounts(27, 0);

  std::vector<Molecule> haloBatch1;
  for (size_t i = 100; i < 115; ++i) {
    haloBatch1.emplace_back(std::array<double, 3>{-0.2, 1.0 + 0.2 * static_cast<double>(i - 100), 1.0},
                            std::array<double, 3>{0.0, 0.0, 0.0}, i);
  }
  autoPas.addHaloParticles(haloBatch1);
  checkCounts(27, 15);

  std::vector<Molecule> haloBatch2;
  for (size_t i = 200; i < 210; ++i) {
    haloBatch2.emplace_back(std::array<double, 3>{10.2, 1.0 + 0.2 * static_cast<double>(i - 200), 1.0},
                            std::array<double, 3>{0.0, 0.0, 0.0}, i);
  }
  // Only even IDs (5 particles) pass predicate.
  autoPas.addHaloParticlesIf(haloBatch2, [](const auto &p) { return p.getID() % 2 == 0; });
  checkCounts(27, 20);
  EXPECT_EQ(countByIterator(autopas::IteratorBehavior::owned | autopas::IteratorBehavior::containerOnly), 27);
  EXPECT_EQ(countByIterator(autopas::IteratorBehavior::halo | autopas::IteratorBehavior::containerOnly), 20);

  // --- Phase 2: Valid neighbor-list state (iteration 1, non-rebuild -> buffers + updateHaloParticle) ---
  LJFunctorType<> functor{autoPas.getCutoff()};
  functor.setParticleProperties(1.0, 1.0);
  autoPas.computeInteractions(&functor);
  auto leaving = autoPas.updateContainer();
  EXPECT_TRUE(leaving.empty());
  // updateContainer marks all halos as dummy while keeping neighbor lists valid.
  checkCounts(27, 0);
  EXPECT_EQ(countByIterator(autopas::IteratorBehavior::owned | autopas::IteratorBehavior::containerOnly), 27);
  EXPECT_EQ(countByIterator(autopas::IteratorBehavior::halo | autopas::IteratorBehavior::containerOnly), 0);

  // Bulk add owned particles into thread-private _particleBuffer.
  std::vector<Molecule> ownedBufferBatch;
  for (size_t i = 300; i < 311; ++i) {
    ownedBufferBatch.emplace_back(std::array<double, 3>{3.0, 3.0, 1.0 + 0.2 * static_cast<double>(i - 300)},
                                  std::array<double, 3>{0.0, 0.0, 0.0}, i);
  }
  autoPas.addParticles(ownedBufferBatch);
  checkCounts(38, 0);
  // Container-only count must remain 27, proving the 11 new particles went into _particleBuffer.
  EXPECT_EQ(countByIterator(autopas::IteratorBehavior::owned | autopas::IteratorBehavior::containerOnly), 27);

  // Bulk add halo particles: haloBatch1 updates existing halo(dummy) cells in place;
  // haloBufferBatch goes into thread-private _haloParticleBuffer.
  autoPas.addHaloParticles(haloBatch1);
  checkCounts(38, 15);
  EXPECT_EQ(countByIterator(autopas::IteratorBehavior::halo | autopas::IteratorBehavior::containerOnly), 15);

  std::vector<Molecule> haloBufferBatch;
  for (size_t i = 400; i < 410; ++i) {
    haloBufferBatch.emplace_back(std::array<double, 3>{10.3, 5.0, 1.0 + 0.2 * static_cast<double>(i - 400)},
                                 std::array<double, 3>{0.0, 0.0, 0.0}, i);
  }
  // Only odd IDs (5 particles) pass predicate.
  autoPas.addHaloParticlesIf(haloBufferBatch, [](const auto &p) { return p.getID() % 2 == 1; });
  checkCounts(38, 20);
  // Container-only halo count must remain 15, proving the 5 new halos went into _haloParticleBuffer.
  EXPECT_EQ(countByIterator(autopas::IteratorBehavior::halo | autopas::IteratorBehavior::containerOnly), 15);
}

/**
 * Verifies that position checks remain enforced when inserting owned or halo particles in valid neighbor-list
 * (buffered) mode.
 */
TEST_F(AutoPasTest, addParticlePositionCheckInBufferModeTest) {
  // Enter iteration 1 (valid neighbor lists -> buffer mode).
  EXPECT_TRUE(autoPas.updateContainer().empty());
  autoPas.addParticle(Molecule({1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, 0));
  LJFunctorType<> functor{autoPas.getCutoff()};
  functor.setParticleProperties(1.0, 1.0);
  autoPas.computeInteractions(&functor);
  EXPECT_TRUE(autoPas.updateContainer().empty());

  // Out-of-box owned particle must still throw in buffer mode.
  EXPECT_ANY_THROW(autoPas.addParticle(Molecule({-0.5, 1.0, 1.0}, {0.0, 0.0, 0.0}, 1)));
  expectedParticles(1, 0);

  // In-box halo particle must still throw in buffer mode.
  EXPECT_ANY_THROW(autoPas.addHaloParticle(Molecule({2.0, 2.0, 2.0}, {0.0, 0.0, 0.0}, 2)));
  expectedParticles(1, 0);
}