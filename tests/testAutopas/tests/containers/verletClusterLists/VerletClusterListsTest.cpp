/**
 * @file VerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "VerletClusterListsTest.h"

#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC06Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterIterationTraversal.h"

namespace VerletClusterListsTest {

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletClusterListsTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  size_t clusterSize = 4;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin, clusterSize);
}

TEST_F(VerletClusterListsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  size_t clusterSize = 4;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, _)).Times(AtLeast(1));
  autopas::VCLClusterIterationTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> verletTraversal(
      &emptyFunctor, clusterSize);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.iteratePairwise(&verletTraversal);
}

TEST_F(VerletClusterListsTest, testAddParticlesAndBuildTwice) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned long numParticles = 271;
  size_t clusterSize = 4;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(
      verletLists, autopas::Particle{}, verletLists.getBoxMin(), verletLists.getBoxMax(), numParticles);

  MockFunctor<Particle> emptyFunctor;
  autopas::VCLClusterIterationTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> verletTraversal(
      &emptyFunctor, clusterSize);
  verletLists.rebuildNeighborLists(&verletTraversal);
  EXPECT_EQ(verletLists.getNumParticles(), numParticles);
  verletLists.rebuildNeighborLists(&verletTraversal);
  EXPECT_EQ(verletLists.getNumParticles(), numParticles);
}

TEST_F(VerletClusterListsTest, testIterator) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  unsigned long numParticles = 271;
  size_t clusterSize = 4;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(
      verletLists, autopas::Particle{}, verletLists.getBoxMin(), verletLists.getBoxMax(), numParticles);

  MockFunctor<Particle> emptyFunctor;
  autopas::VCLClusterIterationTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> verletTraversal(
      &emptyFunctor, clusterSize);
  verletLists.rebuildNeighborLists(&verletTraversal);

  int numParticlesInIterator = 0;
  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    // Make sure that the iterator iterates over real particles, not dummies outside of the domain
    EXPECT_TRUE(autopas::utils::inBox((*it).getR(), min, max));
    numParticlesInIterator++;
  }
  EXPECT_EQ(numParticlesInIterator, numParticles);
}

auto calculateValidPairs(const std::vector<autopas::Particle *> &particles, double cutoffSqr) {
  std::vector<std::pair<autopas::Particle *, autopas::Particle *>> particlePairs;
  for (auto *particlePtr : particles) {
    for (auto *neighborPtr : particles) {
      if (particlePtr != neighborPtr) {
        auto dist = autopas::utils::ArrayMath::sub(particlePtr->getR(), neighborPtr->getR());
        if (autopas::utils::ArrayMath::dot(dist, dist) <= cutoffSqr) {
          particlePairs.emplace_back(particlePtr, neighborPtr);
        }
      }
    }
  }
  return particlePairs;
}

void compareParticlePairs(const std::vector<std::pair<autopas::Particle *, autopas::Particle *>> &first,
                          const std::vector<std::pair<autopas::Particle *, autopas::Particle *>> &second) {
  EXPECT_EQ(first.size(), second.size());
  for (auto pair : first) {
    EXPECT_TRUE(std::find(second.begin(), second.end(), pair) != second.end());
  }
}

TEST_F(VerletClusterListsTest, testNeighborListsValidAfterMovingLessThanHalfSkin) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double cutoffSqr = cutoff * cutoff;
  double skin = 0.2;
  unsigned long numParticles = 271;
  size_t clusterSize = 4;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(
      verletLists, autopas::Particle{}, verletLists.getBoxMin(), verletLists.getBoxMax(), numParticles);
  CollectParticlePairsFunctor functor{cutoff, min, max};
  autopas::VCLClusterIterationTraversal<FPCell, CollectParticlePairsFunctor, autopas::DataLayoutOption::aos, false>
      verletTraversal(&functor, clusterSize);
  verletLists.rebuildNeighborLists(&verletTraversal);

  std::vector<autopas::Particle *> particles;
  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    particles.push_back(&(*it));
  }

  auto referenceParticlePairs = calculateValidPairs(particles, cutoffSqr);

  functor.initTraversal();
  verletLists.iteratePairwise(&verletTraversal);
  functor.endTraversal(false);
  auto calculatedParticlePairs = functor.getParticlePairs();
  compareParticlePairs(referenceParticlePairs, calculatedParticlePairs);

  // Move particles
  int i = -130;
  for (auto it = verletLists.begin(); it.isValid(); ++it, ++i) {
    auto &particle = *it;
    // generate some different directions
    std::array<double, 3> direction = {(double)(i % 2), (double)(i % 3),
                                       (double)(i % numParticles) / (double)numParticles};
    auto offset = autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::normalize(direction), skin / 2.1);
    particle.addR(offset);
    // Upper corner is excluded
    constexpr double smallValue = 0.000001;
    auto x = std::clamp(particle.getR()[0], min[0], max[0] - smallValue);
    auto y = std::clamp(particle.getR()[1], min[1], max[1] - smallValue);
    auto z = std::clamp(particle.getR()[2], min[2], max[2] - smallValue);
    particle.setR({x, y, z});
  }

  auto referenceParticlePairsAfterMove = calculateValidPairs(particles, cutoffSqr);

  functor.initTraversal();
  verletLists.iteratePairwise(&verletTraversal);
  functor.endTraversal(false);
  auto calculatedParticlePairsAfterMove = functor.getParticlePairs();
  compareParticlePairs(referenceParticlePairsAfterMove, calculatedParticlePairsAfterMove);
}

auto getClusterNeighbors(autopas::VerletClusterLists<Particle> &verletLists) {
  std::unordered_map<size_t, std::vector<size_t>> neighbors;
  verletLists.traverseClusters<false>([&neighbors](auto &cluster) {
    auto idFirstParticleInCluster = cluster[0].getID();
    for (const auto &neighborCluster : cluster.getNeighbors()) {
      neighbors[idFirstParticleInCluster].push_back((*neighborCluster)[0].getID());
    }
  });
  return neighbors;
}

size_t getNumNeighbors(const std::unordered_map<size_t, std::vector<size_t>> &neighborList) {
  size_t sum = 0;
  for (auto [id, neighbors] : neighborList) {
    sum += neighbors.size();
  }
  return sum;
}

TEST_F(VerletClusterListsTest, testNewton3NeighborList) {
  std::array<double, 3> min = {0, 0, 0};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.1;
  int numParticles = 2431;
  size_t clusterSize = 4;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(
      verletLists, autopas::Particle{}, verletLists.getBoxMin(), verletLists.getBoxMax(), numParticles);

  MockFunctor<Particle> functor;
  autopas::VCLC06Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> traversalNoN3(&functor,
                                                                                                  clusterSize);
  verletLists.rebuildNeighborLists(&traversalNoN3);
  auto neighborsNoN3 = getClusterNeighbors(verletLists);

  autopas::VCLC06Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> traversalN3(&functor, clusterSize);
  verletLists.rebuildNeighborLists(&traversalN3);
  auto neighborsN3 = getClusterNeighbors(verletLists);

  EXPECT_EQ(getNumNeighbors(neighborsNoN3), getNumNeighbors(neighborsN3) * 2);

  for (auto [idN3, neighbors] : neighborsN3) {
    for (auto idN3Neighbor : neighbors) {
      const auto &idNoN3Neighbors = neighborsNoN3[idN3];
      const auto &idNoN3NeighborNeighbors = neighborsNoN3[idN3Neighbor];
      EXPECT_TRUE(std::find(idNoN3Neighbors.begin(), idNoN3Neighbors.end(), idN3Neighbor) != idNoN3Neighbors.end());
      EXPECT_TRUE(std::find(idNoN3NeighborNeighbors.begin(), idNoN3NeighborNeighbors.end(), idN3) !=
                  idNoN3NeighborNeighbors.end());
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
  size_t clusterSize = 4;
  autopas::VerletClusterLists<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, autopas::Particle{}, min, max,
                                                               numParticles);

  CollectParticlesPerThreadFunctor functor;
  ColoringTraversalWithColorChangeNotify traversal(&functor, clusterSize,
                                                   [&functor](int currentColor) { functor.nextColor(currentColor); });
  functor.initTraversal();
  verletLists.rebuildNeighborLists(&traversal);
  verletLists.iteratePairwise(&traversal);
  functor.endTraversal(true);

  // Check that particles in the same color are only accessed by one thread.
  auto &list = functor._particlesPerThreadPerColor;
  for (int color = 0; color < 8; color++) {
    auto &colorList = list[color];
    for (unsigned long i = 0; i < colorList.size(); i++) {
      for (auto *particlePtr : colorList[i]) {
        for (unsigned long j = i + 1; j < colorList.size(); j++) {
          EXPECT_TRUE(colorList[j].find(particlePtr) == colorList[j].end())
              << particlePtr->toString() << " was accessed by " << i << " and " << j;
        }
      }
    }
  }
}

#endif  // AUTOPAS_OPENMP
} // end namespace VerletClusterListsTest
