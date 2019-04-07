/**
 * @file VerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "VerletClusterClusterTest.h"
#include "autopas/containers/verletClusterLists/VerletClusterCluster.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClusterClusterTraversal.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletClusterClusterTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterCluster<Particle> verletLists(min, max, cutoff, skin);
}

TEST_F(VerletClusterClusterTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterCluster<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, false)).Times(AtLeast(1));
  autopas::VerletClusterClusterTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> dummyTraversal(
      &emptyFunctor);
  verletLists.iteratePairwiseAoS(&emptyFunctor, &dummyTraversal, false);
}
