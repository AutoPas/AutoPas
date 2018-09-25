/**
 * @file VerletListsTest.cpp
 * @author seckler
 * @date 19.04.18
 */

#include "VerletListsTest.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletListsTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);
}

TEST_F(VerletListsTest, testAddParticleNumParticle) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);
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

TEST_F(VerletListsTest, testDeleteAllParticles) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);
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

TEST_F(VerletListsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  auto& list = verletLists.getVerletListsAoS();

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (auto i : list) {
    partners += i.second.size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_F(VerletListsTest, testVerletList) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &mockFunctor);
  verletLists.iteratePairwiseAoS(&mockFunctor, &dummyTraversal, true);

  auto& list = verletLists.getVerletListsAoS();

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (auto i : list) {
    partners += i.second.size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_F(VerletListsTest, testVerletListInSkin) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {1.4, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &mockFunctor);
  verletLists.iteratePairwiseAoS(&mockFunctor, &dummyTraversal, true);

  auto& list = verletLists.getVerletListsAoS();

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (auto i : list) {
    partners += i.second.size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_F(VerletListsTest, testVerletListBuildTwice) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  auto& list = verletLists.getVerletListsAoS();

  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (auto i : list) {
    partners += i.second.size();
  }
  EXPECT_EQ(partners, 1);
}

TEST_F(VerletListsTest, testVerletListBuildFarAway) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {5, 5, 5};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);

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
  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  auto& list = verletLists.getVerletListsAoS();

  ASSERT_EQ(list.size(), 3);
  int partners = 0;
  for (auto i : list) {
    partners += i.second.size();
  }
  ASSERT_EQ(partners, 1);
}

TEST_F(VerletListsTest, testVerletListBuildHalo) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {0.9, 0.9, 0.9};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addHaloParticle(p);
  std::array<double, 3> r2 = {1.1, 1.1, 1.1};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));
  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  auto& list = verletLists.getVerletListsAoS();

  ASSERT_EQ(list.size(), 2);
  EXPECT_EQ(list.size(), 2);
  int partners = 0;
  for (auto i : list) {
    partners += i.second.size();
  }
  ASSERT_EQ(partners, 1);
}

TEST_F(VerletListsTest, testRebuildFrequencyAlways) {
  MockVerletLists<Particle> mockVerletLists({0., 0., 0.}, {10., 10., 10.}, 1., 0.3, 1);

  MockFunctor<Particle, FPCell> emptyFunctor;
  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(4);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
}

TEST_F(VerletListsTest, testRebuildFrequencyEvery3) {
  MockVerletLists<Particle> mockVerletLists({0., 0., 0.}, {10., 10., 10.}, 1., 0.3, 3);

  MockFunctor<Particle, FPCell> emptyFunctor;
  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);

  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);  // 1
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);  // 2
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);  // 3
}

TEST_F(VerletListsTest, testForceRebuild) {
  // generate Velet list with rebuild frequency of 3
  MockVerletLists<Particle> mockVerletLists({0., 0., 0.}, {10., 10., 10.}, 1., 0.3, 3);
  // delegating to parent
  ON_CALL(mockVerletLists, addParticle(_))
      .WillByDefault(Invoke(&mockVerletLists, &MockVerletLists<Particle>::addParticleVerletLists));
  // delegating to parent
  ON_CALL(mockVerletLists, addHaloParticle(_))
      .WillByDefault(Invoke(&mockVerletLists, &MockVerletLists<Particle>::addHaloParticleVerletLists));
  // delegating to parent
  ON_CALL(mockVerletLists, updateContainer())
      .WillByDefault(Invoke(&mockVerletLists, &MockVerletLists<Particle>::updateContainerVerletLists));

  MockFunctor<Particle, FPCell> emptyFunctor;
  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  // check that the second call does not need a rebuild
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  // check that updateContainer() requires a rebuild
  EXPECT_CALL(mockVerletLists, updateContainer());
  mockVerletLists.updateContainer();
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);

  // check that adding particles requires a rebuild
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  EXPECT_CALL(mockVerletLists, addParticle(_));
  mockVerletLists.addParticle(p);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, true);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);

  // check that adding halo particles requires a rebuild
  Particle p2({-0.1, 1.2, 1.1}, {0., 0., 0.}, 2);
  EXPECT_CALL(mockVerletLists, addHaloParticle(_));
  mockVerletLists.addHaloParticle(p2);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal,
                                     true);  // rebuild happens here

  // check that deleting particles requires a rebuild
  /// @todo: reenable once implemented
  /*{
    auto iterator = mockVerletLists.begin();
    iterator.deleteCurrentParticle();
  }
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,&dummyTraversal,
                                     true);  // rebuild happens here
*/
  // check that deleting halo particles requires a rebuild
  mockVerletLists.deleteHaloParticles();
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal,
                                     true);  // rebuild happens here
}

TEST_F(VerletListsTest, testCheckNeighborListsAreValidAfterBuild) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, true)).Times(AtLeast(1));

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  // this will build the verlet list
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal);

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VerletListsTest, testCheckNeighborListsAreValidAfterSmallMove) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  // this will build the verlet list
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.4, 1.1, 1.1});
    }
  }

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VerletListsTest, testCheckNeighborListsAreInvalidAfterMoveLarge) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  // this will build the verlet list
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.6, 1.1, 1.1});
    }
  }

  // check validity - should return true
  EXPECT_FALSE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VerletListsTest, testCheckNeighborListsInvalidMoveFarOutsideCell) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // addtwo particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  // this will build the verlet list
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      // this sets the particle more than skin/2 outside of cell (xmax_cell=2.3)
      iter->setR({2.7, 1.1, 1.1});
    }
  }

  // check validity - should return true
  EXPECT_FALSE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VerletListsTest, testCheckNeighborListsValidMoveLittleOutsideCell) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  MockFunctor<Particle, FPCell> emptyFunctor;

  // add two particles in proper distance
  Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &emptyFunctor);
  // this will build the verlet list
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal);

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
void moveUpdateAndExpectEqual(Container& container, Particle& particle, std::array<double, 3> newPosition) {
  particle.setR(newPosition);
  container.updateHaloParticle(particle);
  {
    auto iter = container.begin();
    auto r = iter->getR();
    EXPECT_EQ(r[0], newPosition[0]);
    EXPECT_EQ(r[1], newPosition[1]);
    EXPECT_EQ(r[2], newPosition[2]);
  }
};

TEST_F(VerletListsTest, testUpdateHaloParticle) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  Particle p({-.1, 10.1, -.1}, {0., 0., 0.}, 1);
  verletLists.addHaloParticle(p);

  // test same position, change velocity
  p.setV({.1, .1, .1});
  verletLists.updateHaloParticle(p);
  {
    auto iter = verletLists.begin();
    auto v = iter->getV();
    for (int i = 0; i < 3; i++) EXPECT_EQ(v[i], 0.1);
  }

  // test different position, same cell
  moveUpdateAndExpectEqual(verletLists, p, {-.05, 10.1, -.1});

  // test different position, neighboring cells
  EXPECT_NO_THROW(moveUpdateAndExpectEqual(verletLists, p, {.05, 10.1, -.1}));
  EXPECT_NO_THROW(moveUpdateAndExpectEqual(verletLists, p, {-.1, 9.95, -.1}));
  EXPECT_NO_THROW(moveUpdateAndExpectEqual(verletLists, p, {-.1, 10.1, .05}));
  EXPECT_NO_THROW(moveUpdateAndExpectEqual(verletLists, p, {-.1, 9.95, .05}));
  EXPECT_NO_THROW(moveUpdateAndExpectEqual(verletLists, p, {.05, 10.1, .05}));
  EXPECT_NO_THROW(moveUpdateAndExpectEqual(verletLists, p, {.05, 9.95, -.1}));
  EXPECT_NO_THROW(moveUpdateAndExpectEqual(verletLists, p, {.05, 9.95, .05}));

  // check for particle with wrong id
  Particle p2({-.1, -.1, -.1}, {0., 0., 0.}, 2);
  EXPECT_ANY_THROW(verletLists.updateHaloParticle(p2));

  // test move far, expect throw
  EXPECT_ANY_THROW(moveUpdateAndExpectEqual(verletLists, p, {3, 3, 3}););

  // test particles at intermediate positions (not at corners)
  Particle p3({-1., 4., 2.}, {0., 0., 0.}, 3);
  verletLists.addHaloParticle(p3);
  EXPECT_NO_THROW(verletLists.updateHaloParticle(p3));
  Particle p4({4., 10.2, 2.}, {0., 0., 0.}, 4);
  verletLists.addHaloParticle(p4);
  EXPECT_NO_THROW(verletLists.updateHaloParticle(p4));
  Particle p5({5., 4., 10.2}, {0., 0., 0.}, 3);
  verletLists.addHaloParticle(p5);
  EXPECT_NO_THROW(verletLists.updateHaloParticle(p5));
}

TEST_F(VerletListsTest, testIsContainerNeeded) {
  std::array<double, 3> boxMin{0, 0, 0};
  std::array<double, 3> boxMax{10, 10, 10};
  double cutoff = 1.;
  double skin = 1.;
  autopas::VerletLists<Particle> container(boxMin, boxMax, cutoff, skin);

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

TEST_F(VerletListsTest, LoadExtractSoA) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  Particle p({-.1, 10.1, -.1}, {0., 0., 0.}, 1);
  verletLists.addHaloParticle(p);

  MockFunctor<Particle, FPCell> mockFunctor;
  autopas::C08Traversal<FPCell, MFunctor, false, true> dummyTraversal({0, 0, 0}, &mockFunctor);

  EXPECT_CALL(mockFunctor, SoALoaderVerlet(_, _, _)).Times(216);  // 6*6*6=216 cells
  EXPECT_CALL(mockFunctor, SoAExtractorVerlet(_, _, _)).Times(216);
  EXPECT_CALL(mockFunctor, SoAFunctor(_, _, _, _, _)).Times(1);
  verletLists.iteratePairwiseSoA(&mockFunctor, &dummyTraversal, true);
}

TEST_F(VerletListsTest, LoadExtractSoALJ) {
  autopas::VerletLists<Particle> verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  Particle p({-.1, 10.1, -.1}, {0., 0., 0.}, 1);
  verletLists.addHaloParticle(p);

  autopas::LJFunctor<Particle, FPCell> ljFunctor;
  autopas::C08Traversal<FPCell, autopas::LJFunctor<Particle, FPCell>, false, true> dummyTraversal({0, 0, 0},
                                                                                                  &ljFunctor);
  verletLists.iteratePairwiseSoA(&ljFunctor, &dummyTraversal, true);
}

TEST_F(VerletListsTest, SoAvsAoSLJ) {
  double cutoff = 2.;
  autopas::VerletLists<Particle> verletLists1({0., 0., 0.}, {10., 10., 10.}, cutoff, 0.3, 3);

  autopas::VerletLists<Particle> verletLists2({0., 0., 0.}, {10., 10., 10.}, cutoff, 0.3, 3);

  RandomGenerator::fillWithParticles(verletLists1, Particle({0., 0., 0.}, {0., 0., 0.}, 0), 100);
  RandomGenerator::fillWithParticles(verletLists2, Particle({0., 0., 0.}, {0., 0., 0.}, 0), 100);

  autopas::LJFunctor<Particle, FPCell> ljFunctor;
  ljFunctor.setGlobals(cutoff, 1, 1, 0);
  autopas::C08Traversal<FPCell, autopas::LJFunctor<Particle, FPCell>, false, true> dummyTraversal({0, 0, 0},
                                                                                                  &ljFunctor);

  verletLists1.iteratePairwiseAoS(&ljFunctor, &dummyTraversal);
  verletLists2.iteratePairwiseSoA(&ljFunctor, &dummyTraversal);

  auto iter1 = verletLists1.begin();
  auto iter2 = verletLists2.begin();

  for (; iter1.isValid() && iter2.isValid(); ++iter1, ++iter2) {
    ASSERT_NEAR(iter1->getR()[0], iter2->getR()[0], fabs(iter1->getR()[0] * 1e-7));
    ASSERT_NEAR(iter1->getR()[1], iter2->getR()[1], fabs(iter1->getR()[1] * 1e-7));
    ASSERT_NEAR(iter1->getR()[2], iter2->getR()[2], fabs(iter1->getR()[2] * 1e-7));

    EXPECT_NEAR(iter1->getF()[0], iter2->getF()[0], fabs(iter1->getF()[0] * 1e-7));
    EXPECT_NEAR(iter1->getF()[1], iter2->getF()[1], fabs(iter1->getF()[1] * 1e-7));
    EXPECT_NEAR(iter1->getF()[2], iter2->getF()[2], fabs(iter1->getF()[2] * 1e-7));
  }
  EXPECT_FALSE(iter1.isValid());
  EXPECT_FALSE(iter2.isValid());
}