/**
 * @file VerletListsTest.cpp
 * @author seckler
 * @date 19.04.18
 */

#include "VerletListsTest.h"

#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/particles/OwnershipState.h"
#include "molecularDynamicsLibrary/LJFunctor.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Invoke;
using ::testing::Values;

TEST_P(VerletListsTest, testVerletListBuildAndIterate) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned int rebuildFrequency = 20;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<ParticleFP64> verletLists(min, max, cutoff, skin, rebuildFrequency,
                                                 autopas::VerletLists<ParticleFP64>::BuildVerletListType::VerletSoA,
                                                 cellSizeFactor);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(1);

  autopas::VLListIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(&emptyFunctor,
                                                                              autopas::DataLayoutOption::aos, true);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);

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
  unsigned int rebuildFrequency = 20;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<ParticleFP64> verletLists(min, max, cutoff, skin, rebuildFrequency,
                                                 autopas::VerletLists<ParticleFP64>::BuildVerletListType::VerletSoA,
                                                 cellSizeFactor);

  std::array<double, 3> r = {1.4, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> mockFunctor;
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::VLListIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(&mockFunctor,
                                                                              autopas::DataLayoutOption::aos, true);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);

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
  unsigned int rebuildFrequency = 20;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<ParticleFP64> verletLists(min, max, cutoff, skin, rebuildFrequency,
                                                 autopas::VerletLists<ParticleFP64>::BuildVerletListType::VerletSoA,
                                                 cellSizeFactor);

  std::array<double, 3> r = {2, 2, 2};
  ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VLListIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(&emptyFunctor,
                                                                              autopas::DataLayoutOption::aos, true);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);
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
  unsigned int rebuildFrequency = 20;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<ParticleFP64> verletLists(min, max, cutoff, skin, rebuildFrequency,
                                                 autopas::VerletLists<ParticleFP64>::BuildVerletListType::VerletSoA,
                                                 cellSizeFactor);

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

  autopas::VLListIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(&emptyFunctor,
                                                                              autopas::DataLayoutOption::aos, true);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);

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
  unsigned int rebuildFrequency = 20;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<ParticleFP64> verletLists(min, max, cutoff, skin, rebuildFrequency,
                                                 autopas::VerletLists<ParticleFP64>::BuildVerletListType::VerletSoA,
                                                 cellSizeFactor);

  std::array<double, 3> r = {0.9, 0.9, 0.9};
  ParticleFP64 p(r, {0., 0., 0.}, 0, autopas::OwnershipState::halo);
  verletLists.addHaloParticle(p);
  std::array<double, 3> r2 = {1.1, 1.1, 1.1};
  ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::VLListIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(&emptyFunctor,
                                                                              autopas::DataLayoutOption::aos, true);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);

  auto &list = verletLists.getVerletListsAoS();

  ASSERT_EQ(list.size(), 2);
  int partners = 0;
  for (const auto &i : list) {
    partners += i.second.size();
  }
  ASSERT_EQ(partners, 1);
}

template <class Container, class Particle_T>
bool moveUpdateAndExpectEqual(Container &container, Particle_T &particle, const std::array<double, 3> &newPosition) {
  particle.setR(newPosition);
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
  autopas::VerletLists<ParticleFP64> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 30,
                                                 autopas::VerletLists<ParticleFP64>::BuildVerletListType::VerletSoA,
                                                 cellSizeFactor);

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

TEST_P(VerletListsTest, LoadExtractSoA) {
  const double cutoff = 2.;
  double skin = 0.3;
  unsigned int rebuildFrequency = 30;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<ParticleFP64> verletLists({0., 0., 0.}, {10., 10., 10.}, cutoff, skin, rebuildFrequency,
                                                 autopas::VerletLists<ParticleFP64>::BuildVerletListType::VerletSoA,
                                                 cellSizeFactor);

  ParticleFP64 p({-.1, 10.1, -.1}, {0., 0., 0.}, 1, autopas::OwnershipState::halo);
  verletLists.addHaloParticle(p);

  MockPairwiseFunctor<ParticleFP64> mockFunctor;

  autopas::VLListIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(&mockFunctor,
                                                                              autopas::DataLayoutOption::soa, false);
  const size_t dimWithHalo = 10 / ((cutoff + skin) * cellSizeFactor) + 2ul;
  const size_t numCells = dimWithHalo * dimWithHalo * dimWithHalo;
  EXPECT_CALL(mockFunctor, SoALoader(testing::An<autopas::FullParticleCell<ParticleFP64> &>(), _, _, _))
      .Times(numCells);
  EXPECT_CALL(mockFunctor, SoAExtractor(testing::An<autopas::FullParticleCell<ParticleFP64> &>(), _, _))
      .Times(numCells);
  EXPECT_CALL(mockFunctor, SoAFunctorVerlet(_, _, _, _)).Times(1);

  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);
}

/**
 * This test is only here to assure stuff compiles with the actual functor and does not crash during execution.
 */
TEST_P(VerletListsTest, LoadExtractSoALJ) {
  const double cutoff = 2.;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Molecule> verletLists(
      {0., 0., 0.}, {10., 10., 10.}, cutoff, 0.3 /*skin*/, 30 /*rebuildFrequency*/,
      autopas::VerletLists<Molecule>::BuildVerletListType::VerletSoA, cellSizeFactor);

  Molecule p({-.1, 10.1, -.1}, {0., 0., 0.}, 1, 0);
  p.setOwnershipState(autopas::OwnershipState::halo);
  verletLists.addHaloParticle(p);
  LJFunctorType<> ljFunctor(cutoff);
  ljFunctor.setParticleProperties(1., 1.);
  autopas::VLListIterationTraversal<FMCell, LJFunctorType<>> verletTraversal(&ljFunctor, autopas::DataLayoutOption::soa,
                                                                             false);

  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);
}

TEST_P(VerletListsTest, SoAvsAoSLJ) {
  const double cutoff = 2.;
  const double cellSizeFactor = GetParam();
  autopas::VerletLists<Molecule> verletLists1({0., 0., 0.}, {10., 10., 10.}, cutoff, 0.01, 30,
                                              autopas::VerletLists<Molecule>::BuildVerletListType::VerletSoA,
                                              cellSizeFactor);

  autopas::VerletLists<Molecule> verletLists2({0., 0., 0.}, {10., 10., 10.}, cutoff, 0.01, 30,
                                              autopas::VerletLists<Molecule>::BuildVerletListType::VerletSoA,
                                              cellSizeFactor);

  Molecule defaultParticle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists1, defaultParticle, verletLists1.getBoxMin(),
                                                                verletLists1.getBoxMax(), 100);
  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists2, defaultParticle, verletLists2.getBoxMin(),
                                                                verletLists2.getBoxMax(), 100);
  LJFunctorType<> ljFunctor(cutoff);
  ljFunctor.setParticleProperties(1., 1.);
  autopas::VLListIterationTraversal<FMCell, LJFunctorType<>> verletTraversal1(&ljFunctor,
                                                                              autopas::DataLayoutOption::aos, false);
  autopas::VLListIterationTraversal<FMCell, LJFunctorType<>> soaTraversal(&ljFunctor, autopas::DataLayoutOption::soa,
                                                                          false);
  verletLists1.rebuildNeighborLists(&verletTraversal1);
  verletLists2.rebuildNeighborLists(&soaTraversal);
  verletLists1.computeInteractions(&verletTraversal1);
  verletLists2.computeInteractions(&soaTraversal);

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