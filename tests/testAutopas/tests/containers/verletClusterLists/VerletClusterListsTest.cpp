/**
 * @file VerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "VerletClusterListsTest.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersColoringTraversal.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversal.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletClusterListsTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin);
}

TEST_F(VerletClusterListsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, false)).Times(AtLeast(1));
  autopas::VerletClustersTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> verletTraversal(
      &emptyFunctor);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);
}

// TODO: Add test that rebuilds the container twice to check if particles are not lost during that

/*
int sumNumClusterNeighbors(
    const std::vector<std::vector<std::vector<autopas::ParticleBase<double, unsigned long> *>>> &neighborList) {
  int sum = 0;
  for (const auto &elem : neighborList) {
    for (const auto &inner : elem) {
      sum += inner.size();
    }
  }
  return sum;
}

TEST_F(VerletClusterListsTest, testVerletListNewton3Build) {
  std::array<double, 3> min = {0, 0, 0};
  std::array<double, 3> max = {8, 8, 8};
  double cutoff = 1.;
  double skin = 0.1;
  autopas::VerletClusterLists<Particle> verletListsNoNewton3(min, max, cutoff, skin);
  autopas::VerletClusterLists<Particle> verletListsNewton3(min, max, cutoff, skin);

  RandomGenerator::fillWithParticles(verletListsNoNewton3, autopas::Particle{}, 100);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = verletListsNoNewton3.begin(); it.isValid(); ++it) {
    verletListsNewton3.addParticle(*it);
  }

  // Generate neighbor list without newton 3.
  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, _)).Times(AtLeast(1));
  autopas::VerletClustersColoringTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> noNewton3Traversal(
      &emptyFunctor);
  verletListsNoNewton3.rebuildNeighborLists(&noNewton3Traversal);
  verletListsNoNewton3.iteratePairwise(&noNewton3Traversal);
  const auto &noNewton3NeighborLists = verletListsNoNewton3.getNeighborLists();
  const auto &noNewton3ClusterIndices = verletListsNoNewton3.getClusterIndexMap();
  int noNewton3NumNeighorClusters = sumNumClusterNeighbors(noNewton3NeighborLists);

  // Generate neighbor list with newton 3.
  autopas::VerletClustersColoringTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> newton3Traversal(
      &emptyFunctor);
  verletListsNewton3.rebuildNeighborLists(&newton3Traversal);
  verletListsNewton3.iteratePairwise(&newton3Traversal);
  const auto &newton3NeighborLists = verletListsNewton3.getNeighborLists();
  const auto &newton3ClusterIndices = verletListsNewton3.getClusterIndexMap();
  int newton3NumNeighorClusters = sumNumClusterNeighbors(newton3NeighborLists);

  // Subtract the neighbors of each cluster with itself.
  unsigned long noNewton3NumWithoutOwn = noNewton3NumNeighorClusters - verletListsNoNewton3.getNumClusters();
  unsigned long newton3NumWithoutOwn = newton3NumNeighorClusters - verletListsNewton3.getNumClusters();

  // Neighbor list without newton 3 should have double the inter-grid connections than with newton 3.
  EXPECT_TRUE(noNewton3NumWithoutOwn % 2 == 0);
  EXPECT_EQ(noNewton3NumWithoutOwn / 2, newton3NumWithoutOwn);

  // Both neighbor lists should contain the same number of grids.
  EXPECT_EQ(noNewton3NeighborLists.size(), newton3NeighborLists.size());

  for (unsigned long gridIndex = 0; gridIndex < noNewton3NeighborLists.size(); gridIndex++) {
    // Every grid from both neighbor lists should contain the same number of clusters.
    EXPECT_EQ(noNewton3NeighborLists[gridIndex].size(), newton3NeighborLists[gridIndex].size());

    for (unsigned long clusterIndex = 0; clusterIndex < noNewton3NeighborLists[gridIndex].size(); clusterIndex++) {
      const auto &noNewton3Neighbors = noNewton3NeighborLists[gridIndex][clusterIndex];
      const auto &newton3Neighbors = newton3NeighborLists[gridIndex][clusterIndex];
      // For each cluster in both neighbor lists, the list with no newton 3 should have more entries.
      EXPECT_GE(noNewton3Neighbors.size(), newton3Neighbors.size());

      for (const auto &neighborPointer : newton3Neighbors) {
        // Each entry of the list with newton 3 should be in the list without newton 3.
        unsigned int neighborIndex = newton3ClusterIndices.at(neighborPointer);
        auto noNewton3It = std::find_if(noNewton3Neighbors.begin(), noNewton3Neighbors.end(),
                                        [neighborIndex, &noNewton3ClusterIndices](auto neighbor) {
                                          return neighborIndex == noNewton3ClusterIndices.at(neighbor);
                                        });
        EXPECT_TRUE(noNewton3It != noNewton3Neighbors.end());
      }
    }
  }
}

#if defined(AUTOPAS_OPENMP)
TEST_F(VerletClusterListsTest, testVerletListColoringTraversalNewton3NoDataRace) {
  std::array<double, 3> min = {0, 0, 0};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.1;
  int numParticles = 5000;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin);

  RandomGenerator::fillWithParticles(verletLists, autopas::Particle{}, numParticles);

  CollectParticlesPerThreadFunctor functor;
  ColoringTraversalWithColorChangeNotify traversal(
      &functor, [](int currentColor) { CollectParticlesPerThreadFunctor::nextColor(currentColor); });
  functor.initTraversal();
  verletLists.rebuildNeighborLists(&traversal);
  verletLists.iteratePairwise(&traversal);
  functor.endTraversal(true);

  // Check that particles in the same color are only accessed by one thread.
  auto &list = functor._particlesPerThreadPerColor;
  for (int color = 0; color < 8; color++) {
    auto &colorList = list[color];
    for (unsigned long i = 0; i < colorList.size(); i++) {
      for (auto particlePtr : colorList[i]) {
        for (unsigned long j = i + 1; j < colorList.size(); j++) {
          EXPECT_TRUE(colorList[j].find(particlePtr) == colorList[j].end())
              << particlePtr->toString() << " was accessed by " << i << " and " << j;
        }
      }
    }
  }
}

#endif  // AUTOPAS_OPENMP
*/