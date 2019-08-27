/**
 * @file VerletListsTest.cpp
 * @author seckler
 * @date 19.04.18
 */

#include "VerletListsTest.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/TraversalVerlet.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Invoke;
using ::testing::Values;

TEST_P(VerletListsTest, testAddParticleNumParticle) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists(
      min, max, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);
  EXPECT_EQ(verletLists.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  EXPECT_EQ(verletLists.getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  EXPECT_EQ(verletLists.getNumParticles(), 2);
}

TEST_P(VerletListsTest, testDeleteAllParticles) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists(
      min, max, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);
  EXPECT_EQ(verletLists.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  EXPECT_EQ(verletLists.getNumParticles(), 2);

  verletLists.deleteAllParticles();
  EXPECT_EQ(verletLists.getNumParticles(), 0);
}

TEST_P(VerletListsTest, testVerletListBuildAndIterate) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists(
      min, max, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(1);

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);

  auto &list = verletLists.getVerletListsAoS();

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (const auto &i : list) {
    partners += i.second.size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_P(VerletListsTest, testVerletListInSkin) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists(
      min, max, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);

  std::array<double, 3> r = {1.4, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> mockFunctor;
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&mockFunctor);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);

  auto &list = verletLists.getVerletListsAoS();

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (const auto &i : list) {
    partners += i.second.size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_P(VerletListsTest, testVerletListBuildTwice) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists(
      min, max, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);
  auto &list = verletLists.getVerletListsAoS();

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (const auto &i : list) {
    partners += i.second.size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_P(VerletListsTest, testVerletListBuildFarAway) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists(
      min, max, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  std::array<double, 3> r3 = {4.5, 4.5, 4.5};
  Particle p3(r3, {0., 0., 0.}, 2);
  verletLists.addParticle(p3);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);

  auto &list = verletLists.getVerletListsAoS();

  ASSERT_EQ(list.size(), 3);
  int partners = 0;
  for (const auto &i : list) {
    partners += i.second.size();
  }
  ASSERT_EQ(partners, 1);
}

TEST_P(VerletListsTest, testVerletListBuildHalo) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists(
      min, max, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);

  std::array<double, 3> r = {0.9, 0.9, 0.9};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addHaloParticle(p);
  std::array<double, 3> r2 = {1.1, 1.1, 1.1};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);

  auto &list = verletLists.getVerletListsAoS();

  ASSERT_EQ(list.size(), 2);
  int partners = 0;
  for (const auto &i : list) {
    partners += i.second.size();
  }
  ASSERT_EQ(partners, 1);
}

TEST_P(VerletListsTest, testCheckNeighborListsAreValidAfterBuild) {
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  // this will build the verlet list
  verletLists.rebuildNeighborLists(&verletTraversal);

  verletLists.iteratePairwise(&verletTraversal);

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_P(VerletListsTest, testCheckNeighborListsAreValidAfterSmallMove) {
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  // this will build the verlet list
  verletLists.rebuildNeighborLists(&verletTraversal);

  verletLists.iteratePairwise(&verletTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.4, 1.1, 1.1});
      break;
    }
  }

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_P(VerletListsTest, testCheckNeighborListsAreInvalidAfterMoveLarge) {
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  // this will build the verlet list
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.6, 1.1, 1.1});
      break;
    }
  }

  // check validity - should return true
  EXPECT_FALSE(verletLists.checkNeighborListsAreValid());
}

TEST_P(VerletListsTest, testCheckNeighborListsInvalidMoveFarOutsideCell) {
  const double cutoff = 2.;
  const double skin = 0.3;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, cutoff, skin,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // add two particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  // this will build the verlet list
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      // this sets the particle more than skin/2 outside of cell (xmax_cell=2.3)
      const double cellLength = 10.0 / (static_cast<int>(10.0 / ((cutoff + skin) * cellSizeFactor)));
      iter->setR({cellLength + skin, 1.1, 1.1});
      break;
    }
  }

  // check validity - should return false
  EXPECT_FALSE(verletLists.checkNeighborListsAreValid());
}

TEST_P(VerletListsTest, testCheckNeighborListsValidMoveLittleOutsideCell) {
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // add two particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> verletTraversal(&emptyFunctor);
  // this will build the verlet list
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      // this sets the particle less than skin/2 outside of cell (xmax_cell=2.3)
      iter->setR({2.4, 1.1, 1.1});
      break;
    }
  }

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

template <class Container, class Particle>
bool moveUpdateAndExpectEqual(Container &container, Particle &particle, std::array<double, 3> newPosition) {
  particle.setR(newPosition);
  /// @todo: uncomment
  bool updated = container.updateHaloParticle(particle);
  if (updated) {
    auto iter = container.begin();
    auto r = iter->getR();
    EXPECT_THAT(r, Eq(newPosition));
  }
  return updated;
}

TEST_P(VerletListsTest, testUpdateHaloParticle) {
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  Particle p({-.1, 10.1, -.1}, {0., 0., 0.}, 1);
  verletLists.addHaloParticle(p);

  // test same position, change velocity
  p.setV({.1, .1, .1});

  EXPECT_TRUE(verletLists.updateHaloParticle(p));

  {
    auto iter = verletLists.begin();
    auto v = iter->getV();
    EXPECT_THAT(v, Each(0.1));
  }

  // test different position, same cell
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {-.05, 10.1, -.1}));

  // test different position, neighboring cells
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {.05, 10.1, -.1}));
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {-.1, 9.95, -.1}));
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {-.1, 10.1, .05}));
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {-.1, 9.95, .05}));
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {.05, 10.1, .05}));
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {.05, 9.95, -.1}));
  EXPECT_TRUE(moveUpdateAndExpectEqual(verletLists, p, {.05, 9.95, .05}));

  // check for particle with wrong id
  Particle p2({-.1, -.1, -.1}, {0., 0., 0.}, 2);
  EXPECT_FALSE(verletLists.updateHaloParticle(p2));

  // test move far, expect throw
  EXPECT_FALSE(moveUpdateAndExpectEqual(verletLists, p, {3, 3, 3}));

  // test particles at intermediate positions (not at corners)
  Particle p3({-1., 4., 2.}, {0., 0., 0.}, 3);
  verletLists.addHaloParticle(p3);
  EXPECT_TRUE(verletLists.updateHaloParticle(p3));
  Particle p4({4., 10.2, 2.}, {0., 0., 0.}, 4);
  verletLists.addHaloParticle(p4);
  EXPECT_TRUE(verletLists.updateHaloParticle(p4));
  Particle p5({5., 4., 10.2}, {0., 0., 0.}, 3);
  verletLists.addHaloParticle(p5);
  EXPECT_TRUE(verletLists.updateHaloParticle(p5));
}

TEST_P(VerletListsTest, testIsContainerUpdateNeeded) {
  std::array<double, 3> boxMin{0, 0, 0};
  std::array<double, 3> boxMax{10, 10, 10};
  double cutoff = 1.;
  double skin = 1.;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> container(
      boxMin, boxMax, cutoff, skin, autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA, cellSizeFactor);

  EXPECT_FALSE(container.isContainerUpdateNeeded());

  Particle p({1, 1, 1}, {0, 0, 0}, 0);
  container.addParticle(p);
  EXPECT_FALSE(container.isContainerUpdateNeeded());

  // Particle moves to different cell -> needs update
  const double cellLength = 10.0 / (static_cast<int>(10.0 / ((cutoff + skin) * cellSizeFactor)));
  container.begin()->setR({1.0 + cellLength, 1.0, 1.0});
  EXPECT_TRUE(container.isContainerUpdateNeeded()) << "cellLength: " << cellLength;

  // Particle moves to halo cell -> needs update
  container.begin()->setR({-1, -1, -1});
  EXPECT_TRUE(container.isContainerUpdateNeeded());
}

TEST_P(VerletListsTest, LoadExtractSoA) {
  const double cutoff = 2.;
  const double skin = 0.3;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, cutoff, skin,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  Particle p({-.1, 10.1, -.1}, {0., 0., 0.}, 1);
  verletLists.addHaloParticle(p);

  MockFunctor<Particle, FPCell> mockFunctor;

  autopas::TraversalVerlet<FPCell, MFunctor, autopas::DataLayoutOption::soa, false> verletTraversal(&mockFunctor);
  const size_t dimWithHalo = 10 / ((cutoff + skin) * cellSizeFactor) + 2ul;
  const size_t numCells = dimWithHalo * dimWithHalo * dimWithHalo;
  EXPECT_CALL(mockFunctor, SoALoaderVerlet(_, _, _)).Times(numCells);
  EXPECT_CALL(mockFunctor, SoAExtractorVerlet(_, _, _)).Times(numCells);
  EXPECT_CALL(mockFunctor, SoAFunctor(_, _, _, _)).Times(1);

  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);
}

/**
 * This test is only here to assure stuff compiles with the actual functor and does not crash during execution.
 */
TEST_P(VerletListsTest, LoadExtractSoALJ) {
  const double cutoff = 2.;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, cutoff, 0.3 /*skin*/,
                                             autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                             cellSizeFactor);

  Particle p({-.1, 10.1, -.1}, {0., 0., 0.}, 1);
  verletLists.addHaloParticle(p);

  autopas::LJFunctor<Particle, FPCell> ljFunctor(cutoff, 1 /*epsilon*/, 1 /*sigma*/, 0 /*shift*/);
  autopas::TraversalVerlet<FPCell, autopas::LJFunctor<Particle, FPCell>, autopas::DataLayoutOption::soa, false>
      verletTraversal(&ljFunctor);

  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);
}

TEST_P(VerletListsTest, SoAvsAoSLJ) {
  const double cutoff = 2.;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Particle> verletLists1({0., 0., 0.}, {10., 10., 10.}, cutoff, 0.3,
                                              autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                              cellSizeFactor);

  autopas::VerletLists<Particle> verletLists2({0., 0., 0.}, {10., 10., 10.}, cutoff, 0.3,
                                              autopas::VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                              cellSizeFactor);

  RandomGenerator::fillWithParticles(verletLists1, Particle({0., 0., 0.}, {0., 0., 0.}, 0), 100);
  RandomGenerator::fillWithParticles(verletLists2, Particle({0., 0., 0.}, {0., 0., 0.}, 0), 100);

  autopas::LJFunctor<Particle, FPCell> ljFunctor(cutoff, 1, 1, 0);
  autopas::TraversalVerlet<FPCell, autopas::LJFunctor<Particle, FPCell>, autopas::DataLayoutOption::aos, false>
      verletTraversal1(&ljFunctor);
  autopas::TraversalVerlet<FPCell, autopas::LJFunctor<Particle, FPCell>, autopas::DataLayoutOption::soa, false>
      soaTraversal(&ljFunctor);
  verletLists1.rebuildNeighborLists(&verletTraversal1);
  verletLists2.rebuildNeighborLists(&soaTraversal);
  verletLists1.iteratePairwise(&verletTraversal1);
  verletLists2.iteratePairwise(&soaTraversal);

  auto iter1 = verletLists1.begin();
  auto iter2 = verletLists2.begin();

  for (; iter1.isValid() && iter2.isValid(); ++iter1, ++iter2) {
    for (unsigned int dim = 0; dim < 3; dim++) {
      ASSERT_NEAR(iter1->getR()[dim], iter2->getR()[dim], fabs(iter1->getR()[dim] * 1e-7));
      EXPECT_NEAR(iter1->getF()[dim], iter2->getF()[dim], fabs(iter1->getF()[dim] * 1e-7));
    }
  }
  EXPECT_FALSE(iter1.isValid());
  EXPECT_FALSE(iter2.isValid());
}

INSTANTIATE_TEST_SUITE_P(Generated, VerletListsTest, Values(1.0, 2.0), VerletListsTest::PrintToStringParamName());