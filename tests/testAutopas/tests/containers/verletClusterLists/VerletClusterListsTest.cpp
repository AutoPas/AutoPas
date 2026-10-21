/**
 * @file VerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "VerletClusterListsTest.h"

#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC06Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterIterationTraversal.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/GenerateValidConfigurations.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletClusterListsTest, VerletListConstructor) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const size_t clusterSize = 4;
  const autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);
}

TEST_F(VerletClusterListsTest, testVerletListBuild) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);

  const std::array<double, 3> r = {2, 2, 2};
  const ParticleFP64 p(r, {0., 0., 0.}, 0);
  verletClusterLists.addParticle(p);
  const std::array<double, 3> r2 = {1.5, 2, 2};
  const ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  verletClusterLists.addParticle(p2);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, _)).Times(AtLeast(1));
  autopas::VCLClusterIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(
      emptyFunctor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletClusterLists.rebuildNeighborLists(&verletTraversal);
  verletClusterLists.computeInteractions(&verletTraversal);
}

TEST_F(VerletClusterListsTest, testAddParticlesAndBuildTwice) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const unsigned long numParticles = 271;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);

  autopas::generators::UniformGenerator::fillWithParticles(
      verletClusterLists, ParticleFP64{}, verletClusterLists.getBoxMin(), verletClusterLists.getBoxMax(), numParticles);

  MPairwiseFunctor emptyFunctor;
  autopas::VCLClusterIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(
      emptyFunctor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletClusterLists.rebuildNeighborLists(&verletTraversal);
  EXPECT_EQ(verletClusterLists.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), numParticles);
  verletClusterLists.rebuildNeighborLists(&verletTraversal);
  EXPECT_EQ(verletClusterLists.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), numParticles);
}

TEST_F(VerletClusterListsTest, testIterator) {
  const std::array<double, 3> boxMin = {1, 1, 1};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.2;
  const unsigned long numParticles = 271;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);

  autopas::generators::UniformGenerator::fillWithParticles(
      verletClusterLists, ParticleFP64{}, verletClusterLists.getBoxMin(), verletClusterLists.getBoxMax(), numParticles);

  MockPairwiseFunctor<ParticleFP64> emptyFunctor;
  autopas::VCLClusterIterationTraversal<FPCell, MPairwiseFunctor> verletTraversal(
      emptyFunctor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletClusterLists.rebuildNeighborLists(&verletTraversal);

  int numParticlesInIterator = 0;
  for (auto it = verletClusterLists.begin(); it.isValid(); ++it) {
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
  autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);

  // Fill the container with random particles and build neighbor lists
  autopas::generators::UniformGenerator::fillWithParticles(
      verletClusterLists, ParticleFP64{}, verletClusterLists.getBoxMin(), verletClusterLists.getBoxMax(), numParticles);
  CollectParticlePairsFunctor functor{cutoff, boxMin, boxMax};
  autopas::VCLClusterIterationTraversal<FPCell, CollectParticlePairsFunctor> verletTraversal(
      functor, clusterSize, autopas::DataLayoutOption::aos, false);
  verletClusterLists.rebuildNeighborLists(&verletTraversal);

  // copy all generated particles
  std::vector<ParticleFP64 *> particles;
  for (auto it = verletClusterLists.begin(); it.isValid(); ++it) {
    particles.push_back(&(*it));
  }

  // manually compute all pairs with a N^2 loop
  const auto referenceParticlePairs = calculateValidPairs(particles, cutoffSqr);

  // Use special functor to collect neighbor interactions from container neighbor lists
  functor.initTraversal();
  verletClusterLists.computeInteractions(&verletTraversal);
  functor.endTraversal(false);
  auto calculatedParticlePairs = functor.getParticlePairs();
  // Check that everything is fine so far
  compareParticlePairs(referenceParticlePairs, calculatedParticlePairs,
                       "After initialization, before particles moved.");

  // Move all particles by "random" values within the box
  int i = -130;
  for (auto it = verletClusterLists.begin(); it.isValid(); ++it, ++i) {
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
  verletClusterLists.computeInteractions(&verletTraversal);
  functor.endTraversal(false);
  const auto calculatedParticlePairsAfterMove = functor.getParticlePairs();
  // final comparison
  compareParticlePairs(referenceParticlePairsAfterMove, calculatedParticlePairsAfterMove, "After random movement");
}

auto getClusterNeighbors(autopas::VerletClusterLists<ParticleFP64> &verletClusterLists) {
  std::unordered_map<size_t, std::vector<size_t>> neighbors;
  verletClusterLists.traverseClusters<false>([&neighbors](auto &cluster) {
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
    autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);

    autopas::generators::UniformGenerator::fillWithParticles(verletClusterLists, ParticleFP64{},
                                                             verletClusterLists.getBoxMin(),
                                                             verletClusterLists.getBoxMax(), numParticles);

    MockPairwiseFunctor<ParticleFP64> functor;
    if (newton3) {
      autopas::VCLC06Traversal<FPCell, MPairwiseFunctor> traversalNoN3(functor, clusterSize,
                                                                       autopas::DataLayoutOption::aos, false);
      verletClusterLists.rebuildNeighborLists(&traversalNoN3);
      neighborsNoN3 = getClusterNeighbors(verletClusterLists);
    } else {
      autopas::VCLC06Traversal<FPCell, MPairwiseFunctor> traversalN3(functor, clusterSize,
                                                                     autopas::DataLayoutOption::aos, true);
      verletClusterLists.rebuildNeighborLists(&traversalN3);
      neighborsN3 = getClusterNeighbors(verletClusterLists);
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
  autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);
  size_t numParticles{0};
  // lower corner of the inner box, thus in the first inner tower
  const ParticleFP64 p0{boxMin, {0., 0., 0.}, numParticles++};
  verletClusterLists.addParticle(p0);
  // just outside the lower corner, thus in the halo
  const ParticleFP64 p1{boxMin - 0.1, {0., 0., 0.}, numParticles++};
  verletClusterLists.addHaloParticle(p1);
  // just inside the upper corner, thus in the last inner tower
  const ParticleFP64 p2{boxMax - 0.1, {0., 0., 0.}, numParticles++};
  verletClusterLists.addParticle(p2);
  // just outside the upper corner, thus in the halo
  const ParticleFP64 p3{boxMax, {0., 0., 0.}, numParticles++};
  verletClusterLists.addHaloParticle(p3);

  // we need 256 to 500 Particles for a 5x5 grid of towers with side length 2
  // so we add as many as we are missing
  for (; numParticles < 257;) {
    const ParticleFP64 pDummy{(boxMax - boxMin) / 2., {0., 0., 0.}, numParticles++};
    verletClusterLists.addParticle(pDummy);
  }

  verletClusterLists.rebuildTowersAndClusters(false);

  const int expectedTowersPerInteractionLength = 2;
  EXPECT_EQ(verletClusterLists.getMaxNumTowersPerInteractionLength(), expectedTowersPerInteractionLength);
  const std::array<double, 2> expectedSideLength = {2., 2.};
  EXPECT_EQ(verletClusterLists.getTowerSideLength(), expectedSideLength);
  const std::array<size_t, 2> expectedTowersPerDim = {5, 5};
  const std::array<size_t, 2> expectedTowersPerDimTotal =
      expectedTowersPerDim + (2ul * expectedTowersPerInteractionLength);
  EXPECT_EQ(verletClusterLists.getTowersPerDimension(), expectedTowersPerDimTotal);

  const std::array<size_t, 2> expectedHaloWidthInTowers{expectedTowersPerInteractionLength,
                                                        expectedTowersPerInteractionLength};
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p0.getR()), expectedHaloWidthInTowers);
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p1.getR()), expectedHaloWidthInTowers - 1ul);
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p2.getR()),
            expectedTowersPerDimTotal - 1ul - expectedHaloWidthInTowers);
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p3.getR()),
            expectedTowersPerDimTotal - expectedHaloWidthInTowers);
}

/**
 * Like testGridAlignment, but in a thin box where the number of towers per interaction length differs between x and y.
 * Each dimension has as many halo tower layers as it needs towers to cover one interaction length.
 */
TEST_F(VerletClusterListsTest, testGridAlignmentDifferentTowersPerInteractionLength) {
  using namespace autopas::utils::ArrayMath::literals;

  const std::array<double, 3> boxMin{0., 0., 0.};
  const std::array<double, 3> boxMax{3., 50., 3.};
  const double cutoff{2.5};
  const double skin{0.25};
  const size_t clusterSize{4};
  autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);
  size_t numParticles{0};
  size_t numOwnedParticles{0};
  // just inside the lower corner, thus in the first inner tower
  const ParticleFP64 p0{boxMin + 0.1, {0., 0., 0.}, numParticles++};
  verletClusterLists.addParticle(p0);
  ++numOwnedParticles;
  // just outside the lower corner, thus in the halo
  ParticleFP64 p1{boxMin - 0.1, {0., 0., 0.}, numParticles++};
  p1.setOwnershipState(autopas::OwnershipState::halo);
  verletClusterLists.addHaloParticle(p1);
  // just inside the upper corner, thus in the last inner tower
  const ParticleFP64 p2{boxMax - 0.1, {0., 0., 0.}, numParticles++};
  verletClusterLists.addParticle(p2);
  ++numOwnedParticles;
  // just outside the upper corner, thus in the halo
  ParticleFP64 p3{boxMax + 0.1, {0., 0., 0.}, numParticles++};
  p3.setOwnershipState(autopas::OwnershipState::halo);
  verletClusterLists.addHaloParticle(p3);

  // A lattice of 4 x 50 x 3 owned particles leads to towers with side lengths [1, 50/35].
  for (size_t x = 0; x < 4; ++x) {
    for (size_t y = 0; y < 50; ++y) {
      for (size_t z = 0; z < 3; ++z) {
        const ParticleFP64 p{{(x + 0.5) * 0.75, y + 0.5, z + 0.5}, {0., 0., 0.}, numParticles++};
        verletClusterLists.addParticle(p);
        ++numOwnedParticles;
      }
    }
  }

  verletClusterLists.rebuildTowersAndClusters(false);

  const auto sideLength = verletClusterLists.getTowerSideLength();
  const std::array<size_t, 2> towersPerInteractionLength{
      static_cast<size_t>(std::ceil((cutoff + skin) / sideLength[0])),
      static_cast<size_t>(std::ceil((cutoff + skin) / sideLength[1]))};
  ASSERT_EQ(towersPerInteractionLength, (std::array<size_t, 2>{3, 2}))
      << "The number of towers per interaction length has to differ between x and y for the test to properly test.";

  const std::array<size_t, 2> expectedTowersPerDim = {3, 35};
  const std::array<size_t, 2> expectedTowersPerDimTotal = expectedTowersPerDim + (towersPerInteractionLength * 2ul);
  EXPECT_EQ(verletClusterLists.getTowersPerDimension(), expectedTowersPerDimTotal);

  const auto &expectedHaloWidthInTowers = towersPerInteractionLength;
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p0.getR()), expectedHaloWidthInTowers);
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p1.getR()), expectedHaloWidthInTowers - 1ul);
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p2.getR()),
            expectedTowersPerDimTotal - 1ul - expectedHaloWidthInTowers);
  EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(p3.getR()),
            expectedTowersPerDimTotal - expectedHaloWidthInTowers);

  size_t numIteratedOwnedParticles{0};
  for (auto iter = verletClusterLists.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    ++numIteratedOwnedParticles;
  }
  EXPECT_EQ(numIteratedOwnedParticles, numOwnedParticles);
}

/**
 * Checks that particles directly at the upper box boundary end up on the correct side of it.
 */
TEST_F(VerletClusterListsTest, testGridAlignmentAtBoxBoundary) {
  using namespace autopas::utils::ArrayMath::literals;

  const double cutoff{1.};
  const double skin{0.2};
  const size_t clusterSize{4};
  // The box length, the number of particles, which sets the number of owned towers per dimension, and that number.
  // With 5.1 / 5, a position at boxMax is rounded into the last owned tower. With 5.4 / 3, the largest position below
  // boxMax is rounded into the first halo tower.
  for (const auto &[boxLength, numParticles, expectedNumOwnedTowers] :
       {std::tuple<double, size_t, size_t>{5.1, 300, 5}, std::tuple<double, size_t, size_t>{5.4, 50, 3}}) {
    SCOPED_TRACE("Box length " + std::to_string(boxLength));
    const std::array<double, 3> boxMin{0., 0., 0.};
    const std::array<double, 3> boxMax{boxLength, boxLength, boxLength};
    autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);
    size_t id{0};
    // directly at the upper corner, thus in the halo
    ParticleFP64 haloAtBoxMax{boxMax, {0., 0., 0.}, id++};
    haloAtBoxMax.setOwnershipState(autopas::OwnershipState::halo);
    verletClusterLists.addHaloParticle(haloAtBoxMax);
    // at the largest position below the upper corner, thus in the last inner tower
    const ParticleFP64 ownedBelowBoxMax{{std::nextafter(boxMax[0], boxMin[0]), std::nextafter(boxMax[1], boxMin[1]),
                                         std::nextafter(boxMax[2], boxMin[2])},
                                        {0., 0., 0.},
                                        id++};
    verletClusterLists.addParticle(ownedBelowBoxMax);
    while (id < numParticles) {
      verletClusterLists.addParticle(ParticleFP64{(boxMax - boxMin) / 2., {0., 0., 0.}, id++});
    }

    verletClusterLists.rebuildTowersAndClusters(false);

    const auto sideLength = verletClusterLists.getTowerSideLength();
    ASSERT_EQ(sideLength,
              (std::array<double, 2>{boxLength / expectedNumOwnedTowers, boxLength / expectedNumOwnedTowers}));
    const auto towersPerInteractionLength = static_cast<size_t>(std::ceil((cutoff + skin) / sideLength[0]));
    const auto towersPerDim = verletClusterLists.getTowersPerDimension();
    EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(haloAtBoxMax.getR()),
              towersPerDim - towersPerInteractionLength);
    EXPECT_EQ(verletClusterLists.getTowerBlock().getTowerIndex2DAtPosition(ownedBelowBoxMax.getR()),
              towersPerDim - towersPerInteractionLength - 1ul);
  }
}

/**
 * Checks that every owned particle interacts with all particles within the cutoff, for all VerletClusterLists
 * traversals, when there are 2 towers per interaction length in the x-dimension and 1 in the y-dimension.
 *
 * In the past, this could lead to "unsafe" mixed halo/owned towers: Halo particles were expected only at either end of
 * the tower. Therefore, a halo particle marked that following clusters only contained halo particles. The interaction
 * length mismatch meant that some towers contained the box boundary in the x/y dimension, meaning that they mixed
 * owned and halo particles more than was expected.
 *
 *            "safe" tower                     "unsafe" tower, before the fix
 *            (halo only in z)                 (y-component covers within y-boundary and outside)
 *                                               y < 50   :  y >= 50
 *         +-------------+                     +----------+----------+
 *         | H  H  H  H  |  halo               |          :          |
 *   z = 3 |.............|                     |..........:..........|
 *         | O  O  O  O  |  owned              | O  O     :          |  [O O D D]  treated as halo
 *         | O  O  O  O  |  owned              |          : H H H H  |  [H H H H]  ends the owned clusters
 *         | O  O  O  O  |  owned              | O        : H H H    |  [O H H H]  owned
 *   z = 0 |.............|                     |..........:..........|
 *         | H  H  H  H  |  halo               |          :          |
 *         +-------------+                     +----------+----------+
 *
 * This led to clusters with owned particles being treated as halo particles.
 *
 * This was fixed in #1224.
 *
 * This test recreates that scenario.
 */
TEST_F(VerletClusterListsTest, testInteractionsDifferentTowersPerInteractionLength) {
  using namespace autopas::utils::ArrayMath::literals;

  const std::array<double, 3> boxMin{0., 0., 0.};
  const std::array<double, 3> boxMax{3., 50., 3.};
  const double cutoff{2.5};
  const double skin{0.25};
  const size_t clusterSize{4};

  const auto fillContainer = [](autopas::VerletClusterLists<ParticleFP64> &verletClusterLists) {
    size_t numParticles{0};
    const auto addOwned = [&](const std::array<double, 3> &pos) {
      verletClusterLists.addParticle(ParticleFP64{pos, {0., 0., 0.}, numParticles++});
    };
    const auto addHalo = [&](const std::array<double, 3> &pos) {
      ParticleFP64 p{pos, {0., 0., 0.}, numParticles++};
      p.setOwnershipState(autopas::OwnershipState::halo);
      verletClusterLists.addHaloParticle(p);
    };
    // Add owned particles to each row of towers in the y-dimension, except for the final row of towers.
    for (int y = 0; y <= 18; ++y) {
      for (const double x : {0.75, 2.25}) {
        addOwned({x, 1. + 2.4 * y, 1.5});
      }
    }

    // For one tower from the last row of towers, we add a one owned particle, 7 halo particles (halo in the
    // y-dimension), then 2 owned particles. The clusters in that tower would then be [owned, halo, halo, halo], [halo,
    // halo, halo, halo], [owned, owned, dummy, dummy]. The last cluster will then get forgotten.
    addOwned({0.75, 49., 0.2});
    for (int z = 0; z < 7; ++z) {
      addHalo({0.75, 50.5, 0.5 + 0.3 * z});
    }
    addOwned({0.75, 49.5, 2.6});
    addOwned({0.4, 48., 2.8});
    // Halo particles below the box in y. Their purpose is only to bring the particle count to 73, which sets the tower
    // side lengths to [1.5, 50/18], i.e. 2 towers per interaction length in x and 1 in y.
    for (int i = 0; i < 25; ++i) {
      addHalo({0.1 + 0.11 * i, -1. - 0.05 * i, 0.1 * i});
    }
  };

  // Only pairs whose first particle is owned, because halo-halo interactions are not calculated.
  const auto ownedPairs = [](std::vector<std::pair<ParticleFP64 *, ParticleFP64 *>> pairs) {
    std::erase_if(pairs, [](const auto &pair) { return not pair.first->isOwned(); });
    return pairs;
  };

  const auto interactionLength = cutoff + skin;
  // All VerletClusterLists configurations, but only AoS as SoA not needed. Only one CSF, as VCL does not used this.
  const auto configurations = generateAllValidConfigurations(
      autopas::InteractionTypeOption::pairwise, {autopas::ContainerOption::verletClusterLists},
      autopas::TraversalOption::getAllOptions(), autopas::LoadEstimatorOption::getAllOptions(),
      {autopas::DataLayoutOption::aos}, autopas::Newton3Option::getAllOptions(), {1.});
  for (const auto &configuration : configurations) {
    SCOPED_TRACE(configuration.toShortString(false));
    autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize,
                                                                 configuration.loadEstimator);
    fillContainer(verletClusterLists);
    CollectParticlePairsFunctor functor{cutoff, boxMin - interactionLength, boxMax + interactionLength};
    const auto traversal =
        autopas::TraversalSelector::generateTraversalFromConfig<ParticleFP64, CollectParticlePairsFunctor>(
            configuration, functor, verletClusterLists.getTraversalSelectorInfo());
    ASSERT_NE(traversal, nullptr) << "The traversal is not applicable.";
    verletClusterLists.rebuildNeighborLists(traversal.get());

    const auto sideLength = verletClusterLists.getTowerSideLength();
    const std::array<size_t, 2> towersPerInteractionLength{
        static_cast<size_t>(std::ceil(interactionLength / sideLength[0])),
        static_cast<size_t>(std::ceil(interactionLength / sideLength[1]))};
    ASSERT_EQ(towersPerInteractionLength, (std::array<size_t, 2>{2, 1}))
        << "The number of towers per interaction length has to differ between x and y for this test to properly test.";

    std::vector<ParticleFP64 *> particles;
    for (auto iter = verletClusterLists.begin(); iter.isValid(); ++iter) {
      particles.push_back(&(*iter));
    }
    const auto referenceParticlePairs = ownedPairs(calculateValidPairs(particles, cutoff * cutoff));

    functor.initTraversal();
    verletClusterLists.computeInteractions(traversal.get());
    functor.endTraversal(configuration.newton3);
    compareParticlePairs(referenceParticlePairs, ownedPairs(functor.getParticlePairs()),
                         configuration.toShortString(false));
  }
}

#if defined(AUTOPAS_USE_OPENMP)
TEST_F(VerletClusterListsTest, testVerletListColoringTraversalNewton3NoDataRace) {
  const std::array<double, 3> boxMin = {0, 0, 0};
  const std::array<double, 3> boxMax = {3, 3, 3};
  const double cutoff = 1.;
  const double skin = 0.1;
  const int numParticles = 5000;
  const size_t clusterSize = 4;
  autopas::VerletClusterLists<ParticleFP64> verletClusterLists(boxMin, boxMax, cutoff, skin, clusterSize);

  autopas::generators::UniformGenerator::fillWithParticles(verletClusterLists, ParticleFP64{}, boxMin, boxMax,
                                                           numParticles);

  CollectParticlesPerThreadFunctor functor;
  ColoringTraversalWithColorChangeNotify traversal(functor, clusterSize,
                                                   [&functor](int currentColor) { functor.nextColor(currentColor); });
  functor.initTraversal();
  verletClusterLists.rebuildNeighborLists(&traversal);
  verletClusterLists.computeInteractions(&traversal);
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