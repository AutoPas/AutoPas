/**
 * @file VerletListsTest.cpp
 * @author seckler
 * @date 19.04.18
 */

#include "VerletListsTest.h"
#include "mocks/MockFunctor.h"
#include "mocks/MockVerletLists.h"

TEST_F(VerletListsTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists(min, max, cutoff, skin);
}

TEST_F(VerletListsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  autopas::Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  autopas::Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;
  verletLists.iteratePairwiseAoS2(&emptyFunctor, true);

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
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  autopas::Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  autopas::Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>>
      mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  verletLists.iteratePairwiseAoS2(&mockFunctor, true);

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
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {1.4, 2, 2};
  autopas::Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {2.5, 2, 2};
  autopas::Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>>
      mockFunctor;
  using ::testing::_;  // anything is ok
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true));

  verletLists.iteratePairwiseAoS2(&mockFunctor, true);

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
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  autopas::Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  autopas::Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;
  verletLists.iteratePairwiseAoS2(&emptyFunctor, true);

  verletLists.iteratePairwiseAoS2(&emptyFunctor, true);

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
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  autopas::Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);

  std::array<double, 3> r2 = {1.5, 2, 2};
  autopas::Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  std::array<double, 3> r3 = {4.5, 4.5, 4.5};
  autopas::Particle p3(r3, {0., 0., 0.}, 2);
  verletLists.addParticle(p3);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;
  verletLists.iteratePairwiseAoS2(&emptyFunctor, true);

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
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {0.9, 0.9, 0.9};
  autopas::Particle p(r, {0., 0., 0.}, 0);
  verletLists.addHaloParticle(p);
  std::array<double, 3> r2 = {1.1, 1.1, 1.1};
  autopas::Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;
  verletLists.iteratePairwiseAoS2(&emptyFunctor, true);

  verletLists.iteratePairwiseAoS2(&emptyFunctor, true);

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
  MockVerletLists<autopas::Particle,
                  autopas::FullParticleCell<autopas::Particle>>
      mockVerletLists({0., 0., 0.}, {10., 10., 10.}, 1., 0.3, 1);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(4);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
}

TEST_F(VerletListsTest, testRebuildFrequencyEvery3) {
  MockVerletLists<autopas::Particle,
                  autopas::FullParticleCell<autopas::Particle>>
      mockVerletLists({0., 0., 0.}, {10., 10., 10.}, 1., 0.3, 3);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);  // 1
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);  // 2
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);  // 3
}

TEST_F(VerletListsTest, testForceRebuild) {
  using ::testing::_;
  using ::testing::Invoke;

  // generate Velet list with rebuild frequency of 3
  MockVerletLists<autopas::Particle,
                  autopas::FullParticleCell<autopas::Particle>>
      mockVerletLists({0., 0., 0.}, {10., 10., 10.}, 1., 0.3, 3);
  // delegating to parent
  ON_CALL(mockVerletLists, addParticle(_))
      .WillByDefault(Invoke(
          &mockVerletLists,
          &MockVerletLists<autopas::Particle,
                           autopas::FullParticleCell<autopas::Particle>>::
              addParticleVerletLists));
  // delegating to parent
  ON_CALL(mockVerletLists, addHaloParticle(_))
      .WillByDefault(Invoke(
          &mockVerletLists,
          &MockVerletLists<autopas::Particle,
                           autopas::FullParticleCell<autopas::Particle>>::
              addHaloParticleVerletLists));
  // delegating to parent
  ON_CALL(mockVerletLists, updateContainer())
      .WillByDefault(Invoke(
          &mockVerletLists,
          &MockVerletLists<autopas::Particle,
                           autopas::FullParticleCell<autopas::Particle>>::
              updateContainerVerletLists));

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;

  // check that the second call does not need a rebuild
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);

  // check that updateContainer() requires a rebuild
  EXPECT_CALL(mockVerletLists, updateContainer());
  mockVerletLists.updateContainer();
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);

  // check that adding particles requires a rebuild
  autopas::Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  EXPECT_CALL(mockVerletLists, addParticle(_));
  mockVerletLists.addParticle(p);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor, true);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,
                                     true);  // rebuild happens here
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(0);

  // check that adding halo particles requires a rebuild
  autopas::Particle p2({-0.1, 1.2, 1.1}, {0., 0., 0.}, 2);
  EXPECT_CALL(mockVerletLists, addHaloParticle(_));
  mockVerletLists.addHaloParticle(p2);
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,
                                     true);  // rebuild happens here

  // check that deleting particles requires a rebuild
  /*{
    auto iterator = mockVerletLists.begin();
    iterator.deleteCurrentParticle();
  }
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,
                                     true);  // rebuild happens here
*/
  // check that deleting halo particles requires a rebuild
  mockVerletLists.deleteHaloParticles();
  EXPECT_CALL(mockVerletLists, updateVerletListsAoS(true)).Times(1);
  mockVerletLists.iteratePairwiseAoS(&emptyFunctor,
                                     true);  // rebuild happens here
}

TEST_F(VerletListsTest, testCheckNeighborListsAreValidAfterBuild) {
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;

  // addtwo particles in proper distance
  autopas::Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  autopas::Particle p2({3.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  // this will build the verlet list
  verletLists.iteratePairwiseAoS2(&emptyFunctor);

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VerletListsTest, testCheckNeighborListsAreValidAfterSmallMove) {
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;

  // addtwo particles in proper distance
  autopas::Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  autopas::Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  // this will build the verlet list
  verletLists.iteratePairwiseAoS2(&emptyFunctor);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.4, 1.1, 1.1});
    }
  }

  // check validity - should return true
  EXPECT_TRUE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VerletListsTest, testCheckNeighborListsAreInvalidAfterMoveLarge) {
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;

  // addtwo particles in proper distance
  autopas::Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  autopas::Particle p2({3.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  // this will build the verlet list
  verletLists.iteratePairwiseAoS2(&emptyFunctor);

  for (auto iter = verletLists.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      iter->setR({1.6, 1.1, 1.1});
    }
  }

  // check validity - should return true
  EXPECT_FALSE(verletLists.checkNeighborListsAreValid());
}

TEST_F(VerletListsTest, testCheckNeighborListsInvalidMoveFarOutsideCell) {
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;

  // addtwo particles in proper distance
  autopas::Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  autopas::Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  // this will build the verlet list
  verletLists.iteratePairwiseAoS2(&emptyFunctor);

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
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  autopas::Functor<autopas::Particle,
                   autopas::FullParticleCell<autopas::Particle>>
      emptyFunctor;

  // add two particles in proper distance
  autopas::Particle p({1.1, 1.1, 1.1}, {0., 0., 0.}, 1);
  verletLists.addParticle(p);
  autopas::Particle p2({7.5, 1.1, 1.1}, {0., 0., 0.}, 2);
  verletLists.addParticle(p2);

  // this will build the verlet list
  verletLists.iteratePairwiseAoS2(&emptyFunctor);

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
void moveUpdateAndExpectEqual(Container& container, Particle& particle,
                              std::array<double, 3> newPosition) {
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
  autopas::VerletLists<autopas::Particle,
                       autopas::FullParticleCell<autopas::Particle>>
      verletLists({0., 0., 0.}, {10., 10., 10.}, 2., 0.3, 3);

  autopas::Particle p({-.1, 10.1, -.1}, {0., 0., 0.}, 1);
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
  autopas::Particle p2({-.1, -.1, -.1}, {0., 0., 0.}, 2);
  EXPECT_ANY_THROW(verletLists.updateHaloParticle(p2));

  // test move far, expect throw
  EXPECT_ANY_THROW(moveUpdateAndExpectEqual(verletLists, p, {3, 3, 3}););
}