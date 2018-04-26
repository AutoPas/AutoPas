/**
 * @file VerletListsTest.cpp
 * @author seckler
 * @date 19.04.18
 */

#include "VerletListsTest.h"
#include "mocks/MockFunctor.h"

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
