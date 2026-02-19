/**
 * @file VerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "VerletClusterListsTest.h"

#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC06Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC08Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterIterationTraversal.h"
#include "autopas/utils/WrapOpenMP.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletClusterListsTest, VerletListConstructor) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const size_t clusterSize = 4;
  const autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);
}

TEST_F(VerletClusterListsTest, testVerletListBuild) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);

  const std::array<double, 3> r = {2, 2, 2};
  const ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  const std::array<double, 3> r2 = {1.5, 2, 2};
  const ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, _)).Times(AtLeast(1));
  autopas::VCLClusterIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(
      &emptyFunctor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletLists.rebuildNeighborLists(&verletTraversal);
  verletLists.computeInteractions(&verletTraversal);
}

TEST_F(VerletClusterListsTest, testAddParticlesAndBuildTwice) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const unsigned long numParticles = 271;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);

  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists, ParticleFP64{}, verletLists.getBoxMin(),
                                                                verletLists.getBoxMax(), numParticles);

  MPairwiseFunctor emptyFunctor;
  autopas::VCLClusterIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(
      &emptyFunctor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletLists.rebuildNeighborLists(&verletTraversal);
  EXPECT_EQ(verletLists.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), numParticles);
  verletLists.rebuildNeighborLists(&verletTraversal);
  EXPECT_EQ(verletLists.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), numParticles);
}

TEST_F(VerletClusterListsTest, testIterator) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const unsigned long numParticles = 271;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);

  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists, ParticleFP64{}, verletLists.getBoxMin(),
                                                                verletLists.getBoxMax(), numParticles);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  autopas::VCLClusterIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(
      &emptyFunctor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletLists.rebuildNeighborLists(&verletTraversal);

  int numParticlesInIterator = 0;
  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    // Make sure that the iterator iterates over real particles, not dummies outside the domain
    EXPECT_TRUE(autopas::utils::inBox((*it).getR(), boxMin, boxMax));
    numParticlesInIterator++;
  }
  EXPECT_EQ(numParticlesInIterator, numParticles);
}

auto calculateValidPairs(const std::vector<ParticleFP64 *> &particles, double cutoffSqr) {
  using namespace autopas::utils::ArrayMath::literals;

  std::vector<std::pair<ParticleFP64 *, ParticleFP64 *>> particlePairs;
  particlePairs.reserve(particlePairs.size() * (particlePairs.size() / 2));
  for (auto *particlePtr : particles) {
    for (auto *neighborPtr : particles) {
      if (particlePtr != neighborPtr) {
        auto dist = particlePtr->getR() - neighborPtr->getR();
        if (autopas::utils::ArrayMath::dot(dist, dist) <= cutoffSqr) {
          particlePairs.emplace_back(particlePtr, neighborPtr);
        }
      }
    }
  }
  return particlePairs;
}

void compareParticlePairs(const std::vector<std::pair<ParticleFP64 *, ParticleFP64 *>> &expected,
                          const std::vector<std::pair<ParticleFP64 *, ParticleFP64 *>> &actual,
                          const std::string &errMsg) {
  using autopas::utils::ArrayUtils::operator<<;
  EXPECT_EQ(expected.size(), actual.size()) << errMsg << "\n"
                                            << "Number of pairwise interactions differs.";
  EXPECT_THAT(actual, ::testing::UnorderedElementsAreArray(expected)) << errMsg;
}

/**
 * Fill a VCL container with randomly placed particles and compare pairwise interactions before and after
 * random movement to N^2 calculated interactions.
 */
TEST_F(VerletClusterListsTest, testNeighborListsValidAfterMovingLessThanHalfSkin) {
  using namespace autopas::utils::ArrayMath::literals;

  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double cutoffSqr = cutoff * cutoff;
  const double skin = 0.2;
  const unsigned long numParticles = 30;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);

  // Fill the container with random particles and build neighbor lists
  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists, ParticleFP64{}, verletLists.getBoxMin(),
                                                                verletLists.getBoxMax(), numParticles);
  CollectParticlePairsFunctor functor{cutoff, boxMin, boxMax};
  autopas::VCLClusterIterationTraversal<FPCell, CollectParticlePairsFunctor> verletTraversal(
      &functor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletLists.rebuildNeighborLists(&verletTraversal);

  // copy all generated particles
  std::vector<ParticleFP64 *> particles;
  for (auto it = verletLists.begin(); it.isValid(); ++it) {
    particles.push_back(&(*it));
  }

  // manually compute all pairs with a N^2 loop
  const auto referenceParticlePairs = calculateValidPairs(particles, cutoffSqr);

  // Use special functor to collect neighbor interactions from container neighbor lists
  functor.initTraversal();
  verletLists.computeInteractions(&verletTraversal);
  functor.endTraversal(false);
  auto calculatedParticlePairs = functor.getParticlePairs();
  // Check that everything is fine so far
  compareParticlePairs(referenceParticlePairs, calculatedParticlePairs,
                       "After initialization, before particles moved.");

  // Move all particles by "random" values within the box
  int i = -130;
  for (auto it = verletLists.begin(); it.isValid(); ++it, ++i) {
    auto &particle = *it;
    // generate some different directions
    const std::array<double, 3> direction = {static_cast<double>(i % 2), static_cast<double>(i % 3),
                                             static_cast<double>(i % numParticles) / static_cast<double>(numParticles)};
    const auto offset = autopas::utils::ArrayMath::normalize(direction) * (skin / 2.1);
    particle.addR(offset);
    // Clamp the new position to the domain boundaries
    // Upper corner is excluded
    constexpr double smallValue = 0.000001;
    const auto x = std::clamp(particle.getR()[0], boxMin[0], boxMax[0] - smallValue);
    const auto y = std::clamp(particle.getR()[1], boxMin[1], boxMax[1] - smallValue);
    const auto z = std::clamp(particle.getR()[2], boxMin[2], boxMax[2] - smallValue);
    particle.setR({x, y, z});
  }

  // generate a new reference
  const auto referenceParticlePairsAfterMove = calculateValidPairs(particles, cutoffSqr);

  // generate new pairwise interactions
  functor.initTraversal();
  verletLists.computeInteractions(&verletTraversal);
  functor.endTraversal(false);
  const auto calculatedParticlePairsAfterMove = functor.getParticlePairs();
  // final comparison
  compareParticlePairs(referenceParticlePairsAfterMove, calculatedParticlePairsAfterMove, "After random movement");
}

auto getClusterNeighbors(autopas::VerletClusterLists<ParticleFP64> &verletLists) {
  std::unordered_map<size_t, std::vector<size_t>> neighbors;
  verletLists.traverseClusters<false>([&neighbors](auto &cluster) {
    auto idFirstParticleInCluster = cluster[0].getID();
    for (const auto &neighborCluster : *cluster.getNeighbors()) {
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
  const std::array<double, 3> boxMin = {0, 0, 0};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.1;
  const int numParticles = 2431;
  const size_t clusterSize = 4;

  std::unordered_map<size_t, std::vector<size_t>> neighborsNoN3{}, neighborsN3{};

  for (bool newton3 : {true, false}) {
    autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);

    autopasTools::generators::UniformGenerator::fillWithParticles(verletLists, ParticleFP64{}, verletLists.getBoxMin(),
                                                                  verletLists.getBoxMax(), numParticles);

    MockPairwiseFunctor<ParticleFP64> functor;
    if (newton3) {
      autopas::VCLC06Traversal<FPCell, MPairwiseFunctor> traversalNoN3(&functor, clusterSize,
                                                                       autopas::DataLayoutOption::aos, false);
      verletLists.rebuildNeighborLists(&traversalNoN3);
      neighborsNoN3 = getClusterNeighbors(verletLists);
    } else {
      autopas::VCLC06Traversal<FPCell, MPairwiseFunctor> traversalN3(&functor, clusterSize,
                                                                     autopas::DataLayoutOption::aos, true);
      verletLists.rebuildNeighborLists(&traversalN3);
      neighborsN3 = getClusterNeighbors(verletLists);
    }
  }

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

TEST_F(VerletClusterListsTest, testGridAlignment) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::ArrayUtils::static_cast_copy_array;

  const std::array<double, 3> boxMin{0., 0., 0.};
  const std::array<double, 3> boxMax{10., 10., 10.};
  const double cutoff{2.};
  const double skin{0.05};
  const size_t clusterSize{4};
  autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);
  size_t numParticles{0};
  // lower corner of the inner box, thus in the first inner tower
  const ParticleFP64 p0{boxMin, {0., 0., 0.}, numParticles++};
  verletLists.addParticle(p0);
  // just outside the lower corner, thus in the halo
  const ParticleFP64 p1{boxMin - 0.1, {0., 0., 0.}, numParticles++};
  verletLists.addHaloParticle(p1);
  // just inside the upper corner, thus in the last inner tower
  const ParticleFP64 p2{boxMax - 0.1, {0., 0., 0.}, numParticles++};
  verletLists.addParticle(p2);
  // just outside the upper corner, thus in the halo
  const ParticleFP64 p3{boxMax, {0., 0., 0.}, numParticles++};
  verletLists.addHaloParticle(p3);

  // we need 256 to 500 Particles for a 5x5 grid of towers with side length 2
  // so we add as many as we are missing
  for (; numParticles < 257;) {
    const ParticleFP64 pDummy{(boxMax - boxMin) / 2., {0., 0., 0.}, numParticles++};
    verletLists.addParticle(pDummy);
  }

  verletLists.rebuildTowersAndClusters(false);

  const int expectedTowersPerInteractionLength = 2;
  EXPECT_EQ(verletLists.getNumTowersPerInteractionLength(), expectedTowersPerInteractionLength);
  const std::array<double, 2> expectedSideLength = {2., 2.};
  EXPECT_EQ(verletLists.getTowerSideLength(), expectedSideLength);
  const std::array<size_t, 2> expectedTowersPerDim = {5, 5};
  const std::array<size_t, 2> expectedTowersPerDimTotal =
      expectedTowersPerDim + (2ul * expectedTowersPerInteractionLength);
  EXPECT_EQ(verletLists.getTowersPerDimension(), expectedTowersPerDimTotal);

  const std::array<size_t, 2> expectedHaloWidthInTowers{expectedTowersPerInteractionLength,
                                                        expectedTowersPerInteractionLength};
  EXPECT_EQ(verletLists.getTowerBlock().getTowerIndex2DAtPosition(p0.getR()), expectedHaloWidthInTowers);
  EXPECT_EQ(verletLists.getTowerBlock().getTowerIndex2DAtPosition(p1.getR()), expectedHaloWidthInTowers - 1ul);
  EXPECT_EQ(verletLists.getTowerBlock().getTowerIndex2DAtPosition(p2.getR()),
            expectedTowersPerDimTotal - 1ul - expectedHaloWidthInTowers);
  EXPECT_EQ(verletLists.getTowerBlock().getTowerIndex2DAtPosition(p3.getR()),
            expectedTowersPerDimTotal - expectedHaloWidthInTowers);
}

#if defined(AUTOPAS_USE_OPENMP)
TEST_F(VerletClusterListsTest, testVerletListColoringTraversalNewton3NoDataRace) {
  const std::array<double, 3> boxMin = {0, 0, 0};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.1;
  const int numParticles = 5000;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletLists(boxMin, boxMax, cutoff, skin, clusterSize);

  autopasTools::generators::UniformGenerator::fillWithParticles(verletLists, ParticleFP64{}, boxMin, boxMax,
                                                                numParticles);

  CollectParticlesPerThreadFunctor functor;
  ColoringTraversalWithColorChangeNotify traversal(&functor, clusterSize,
                                                   [&functor](int currentColor) { functor.nextColor(currentColor); });
  functor.initTraversal();
  verletLists.rebuildNeighborLists(&traversal);
  verletLists.computeInteractions(&traversal);
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

#endif  // AUTOPAS_USE_OPENMP