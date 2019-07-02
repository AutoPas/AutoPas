/**
 * @file VarVerletListsTest.cpp
 * @author
 * @date 21.05.19
 *
 * Mostly copied from VerletListsTest.cpp
 */

#include "VarVerletListsTest.h"
#include "autopas/containers/verletListsCellBased/verletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/neighborLists/VerletNeighborListAsBuild.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VarVerletTraversalAsBuild.h"
#include "autopas/options/DataLayoutOption.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Invoke;
using ::testing::Values;

TEST_F(VarVerletListsTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);
}

TEST_F(VarVerletListsTest, testAddParticleNumParticle) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);
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

TEST_F(VarVerletListsTest, testDeleteAllParticles) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);
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

TEST_F(VarVerletListsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletList) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&mockFunctor);
  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.iteratePairwise(&mockFunctor, &dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListInSkin) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {1.4, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&mockFunctor);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.iteratePairwise(&mockFunctor, &dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListBuildTwice) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal);

  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListBuildFarAway) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);

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
  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);
  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListBuildHalo) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {0.9, 0.9, 0.9};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addHaloParticle(p);
  std::array<double, 3> r2 = {1.1, 1.1, 1.1};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal);

  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testCheckNeighborListsAreValidAfterBuild) {
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 1);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  // this will build the verlet list
  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal);

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VarVerletListsTest, testCheckNeighborListsAreValidAfterSmallMove) {
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 1);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  // this will build the verlet list
  verletLists.rebuildNeighborLists(&dummyTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.4, 1.1, 1.1});
    }
  }

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VarVerletListsTest, testCheckNeighborListsAreInvalidAfterMoveLarge) {
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 1);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  // this will build the verlet list
  verletLists.rebuildNeighborLists(&dummyTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.6, 1.1, 1.1});
    }
  }

  // check validity - should return false
  EXPECT_FALSE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VarVerletListsTest, testCheckNeighborListsInvalidMoveFarOutsideCell) {
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 1);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  // this will build the verlet list
  verletLists.rebuildNeighborLists(&dummyTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      // this sets the particle more than skin/2 outside of cell (xmax_cell=2.3)
      iter->setR({2.7, 1.1, 1.1});
    }
  }

  // check validity - should return true
  EXPECT_FALSE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VarVerletListsTest, testCheckNeighborListsValidMoveLittleOutsideCell) {
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 1);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // add two particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::VarVerletTraversalAsBuild<FPCell, autopas::Particle, MFunctor, autopas::DataLayoutOption::aos, true>
      dummyTraversal(&emptyFunctor);

  // this will build the verlet list
  verletLists.rebuildNeighborLists(&dummyTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      // this sets the particle less than skin/2 outside of cell (xmax_cell=2.3)
      iter->setR({2.4, 1.1, 1.1});
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

TEST_F(VarVerletListsTest, testUpdateHaloParticle) {
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 1);

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

TEST_F(VarVerletListsTest, testIsContainerNeeded) {
  std::array<double, 3> boxMin{0, 0, 0};
  std::array<double, 3> boxMax{10, 10, 10};
  double cutoff = 1.;
  double skin = 1.;
  autopas::VarVerletLists<Particle, autopas::VerletNeighborListAsBuild<Particle>> container(boxMin, boxMax, cutoff,
                                                                                            skin);

  EXPECT_FALSE(container.isContainerUpdateNeeded());

  Particle p({1, 1, 1}, {0, 0, 0}, 0);
  container.addParticle(p);
  EXPECT_FALSE(container.isContainerUpdateNeeded());

  // Particle moves to different cell -> needs update
  container.begin()->setR({2.5, 1, 1});
  EXPECT_TRUE(container.isContainerUpdateNeeded());

  // Particle moves to halo cell -> needs update
  container.begin()->setR({-1, -1, -1});
  EXPECT_TRUE(container.isContainerUpdateNeeded());
}