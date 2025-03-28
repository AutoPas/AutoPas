/**
 * @file VarVerletListsTest.cpp
 * @author
 * @date 21.05.19
 *
 * Mostly copied from VerletListsTest.cpp
 */

#include "VarVerletListsTest.h"

#include "autopas/containers/verletListsCellBased/varVerletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/traversals/VVLAsBuildTraversal.h"
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
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);
}

TEST_F(VarVerletListsTest, testAddParticleNumParticle) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);
  EXPECT_EQ(verletLists.size(), 0);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  EXPECT_EQ(verletLists.size(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  EXPECT_EQ(verletLists.size(), 2);
}

TEST_F(VarVerletListsTest, testDeleteAllParticles) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);
  EXPECT_EQ(verletLists.size(), 0);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);
  EXPECT_EQ(verletLists.size(), 2);

  verletLists.deleteAllParticles();
  EXPECT_EQ(verletLists.size(), 0);
}

TEST_F(VarVerletListsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VVLAsBuildTraversal<FPCell, ParticleFP64, MPairwiseFunctor> dummyTraversal(
      &emptyFunctor, autopas::DataLayoutOption::aos, true);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.computeInteractions(&dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletList) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::VVLAsBuildTraversal<FPCell, ParticleFP64, MPairwiseFunctor> dummyTraversal(
      &mockFunctor, autopas::DataLayoutOption::aos, true);
  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.computeInteractions(&dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListInSkin) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);

  std::array<double, 3> r = {1.4, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::VVLAsBuildTraversal<FPCell, ParticleFP64, MPairwiseFunctor> dummyTraversal(
      &mockFunctor, autopas::DataLayoutOption::aos, true);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.computeInteractions(&dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListBuildTwice) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VVLAsBuildTraversal<FPCell, ParticleFP64, MPairwiseFunctor> dummyTraversal(
      &emptyFunctor, autopas::DataLayoutOption::aos, true);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.computeInteractions(&dummyTraversal);

  verletLists.computeInteractions(&dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListBuildFarAway) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  std::array<double, 3> r3 = {4.5, 4.5, 4.5};
  ParticleFP64 p3(r3, {0., 0., 0.}, 2);
  verletLists.addParticle(p3);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));
  autopas::VVLAsBuildTraversal<FPCell, ParticleFP64, MPairwiseFunctor> dummyTraversal(
      &emptyFunctor, autopas::DataLayoutOption::aos, true);
  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.computeInteractions(&dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

TEST_F(VarVerletListsTest, testVerletListBuildHalo) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      min, max, cutoff, skin, rebuildFrequency);

  std::array<double, 3> r = {0.9, 0.9, 0.9};
  ParticleFP64 p(r, {0., 0., 0.}, 0, autopas::OwnershipState::halo);
  verletLists.addHaloParticle(p);
  std::array<double, 3> r2 = {1.1, 1.1, 1.1};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VVLAsBuildTraversal<FPCell, ParticleFP64, MPairwiseFunctor> dummyTraversal(
      &emptyFunctor, autopas::DataLayoutOption::aos, true);

  verletLists.rebuildNeighborLists(&dummyTraversal);
  verletLists.computeInteractions(&dummyTraversal);

  verletLists.computeInteractions(&dummyTraversal);

  EXPECT_EQ(verletLists.getNumberOfNeighborPairs(), 1);
}

template <class Container, class ParticleFP64>
bool moveUpdateAndExpectEqual(Container &container, ParticleFP64 &particle, std::array<double, 3> newPosition) {
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
  autopas::VarVerletLists<ParticleFP64, autopas::VerletNeighborListAsBuild<ParticleFP64>> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 1);

  ParticleFP64 p({-.1, 10.1, -.1}, {0., 0., 0.}, 1, autopas::OwnershipState::halo);
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
  ParticleFP64 p2({-.1, -.1, -.1}, {0., 0., 0.}, 2, autopas::OwnershipState::halo);
  EXPECT_FALSE(verletLists.updateHaloParticle(p2));

  // test move far, expect throw
  EXPECT_FALSE(moveUpdateAndExpectEqual(verletLists, p, {3, 3, 3}));

  // test particles at intermediate positions (not at corners)
  ParticleFP64 p3({-1., 4., 2.}, {0., 0., 0.}, 3, autopas::OwnershipState::halo);
  verletLists.addHaloParticle(p3);
  EXPECT_TRUE(verletLists.updateHaloParticle(p3));
  ParticleFP64 p4({4., 10.2, 2.}, {0., 0., 0.}, 4, autopas::OwnershipState::halo);
  verletLists.addHaloParticle(p4);
  EXPECT_TRUE(verletLists.updateHaloParticle(p4));
  ParticleFP64 p5({5., 4., 10.2}, {0., 0., 0.}, 3, autopas::OwnershipState::halo);
  verletLists.addHaloParticle(p5);
  EXPECT_TRUE(verletLists.updateHaloParticle(p5));
}
