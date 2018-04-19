/**
 * @file VerletListsTest.cpp
 * @author seckler
 * @date 19.04.18
 */

#include "VerletListsTest.h"

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

  auto& list = verletLists.getVerletLists();

  EXPECT_EQ(list.size(), 2);
  ASSERT_EQ(list[0].size() + list[1].size(), 1);
  if (list[0].size() == 1) {
    EXPECT_EQ(list[0][0], 1);
  } else if (list[1].size() == 1) {
    EXPECT_EQ(list[1][0], 0);
  }
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

  auto& list = verletLists.getVerletLists();

  EXPECT_EQ(list.size(), 2);
  ASSERT_EQ(list[0].size() + list[1].size(), 1);
  if (list[0].size() == 1) {
    EXPECT_EQ(list[0][0], 1);
  } else if (list[1].size() == 1) {
    EXPECT_EQ(list[1][0], 0);
  }
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

  auto& list = verletLists.getVerletLists();

  EXPECT_EQ(list.size(), 3);
  EXPECT_EQ(list[0].size() + list[1].size() + list[2].size(), 1);

}
