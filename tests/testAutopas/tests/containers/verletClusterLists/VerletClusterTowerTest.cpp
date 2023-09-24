/**
 * @file VerletClusterTowerTest.cpp
 * @author humig
 * @date 05.02.2020
 */

#include "VerletClusterTowerTest.h"

#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "testingHelpers/commonTypedefs.h"

template <class Particle>
using ClusterTower = autopas::internal::ClusterTower<Particle>;

TEST_F(VerletClusterTowerTest, testAddParticle) {
  ClusterTower<Particle> tower(4);
  Particle p1;
  Particle p2;

  tower.addParticle(p1);
  tower.addParticle(p2);

  EXPECT_EQ(tower.getNumActualParticles(), 2);
}

static void testClusterGenerationAndDummies(size_t clusterSize) {
  // Test for different number of particles
  for (size_t numParticles : {29, 64, 1, 2, 4, 0}) {
    ClusterTower<Particle> tower(clusterSize);

    for (size_t i = 0; i < numParticles; i++) {
      tower.addParticle(Particle{{0.0, 0.0, (double)i}, {0, 0, 0}, i});
    }
    EXPECT_EQ(tower.getNumActualParticles(), numParticles);

    tower.generateClusters();
    // Check if number of actual particles is still the same, so dummies are not counted as real particles
    EXPECT_EQ(tower.getNumActualParticles(), numParticles);

    auto numParticlesLastCluster = numParticles % clusterSize;
    if (numParticlesLastCluster == 0) numParticlesLastCluster = clusterSize;
    EXPECT_EQ(tower.getNumTailDummyParticles(), clusterSize - numParticlesLastCluster);

    EXPECT_EQ(tower.getNumClusters(), (numParticles + tower.getNumTailDummyParticles()) / clusterSize);

    // Check if the right particles are within each cluster
    for (size_t clusterIndex = 0; clusterIndex < tower.getNumClusters(); clusterIndex++) {
      const auto &cluster = tower.getCluster(clusterIndex);
      for (size_t particleIndex = 0; particleIndex < clusterSize; particleIndex++) {
        // Since the last particle is copied in the last cluster to fill it up, use minimum of expected index and
        // maximal index.
        unsigned long expectedParticleID = std::min(clusterIndex * clusterSize + particleIndex, numParticles - 1);
        EXPECT_EQ(expectedParticleID, cluster[particleIndex].getID());
      }
    }

    // Skip generating neighbor lists

    // Check if dummy particles are filled in correctly. (Dummy particles always have ID
    // std::numeric_limits<size_t>::max(), filled up particles have ID>0.
    tower.setDummyValues(0, 0);
    const auto &lastCluster = tower.getCluster(tower.getNumClusters() - 1);
    for (size_t i = 1; i <= tower.getNumTailDummyParticles(); i++) {
      EXPECT_EQ(lastCluster[clusterSize - i].getID(), std::numeric_limits<size_t>::max());
    }
  }
}

TEST_F(VerletClusterTowerTest, testClusterGenerationAndDummies) {
  // Test with different cluster sizes
  testClusterGenerationAndDummies(1);
  testClusterGenerationAndDummies(2);
  testClusterGenerationAndDummies(4);
}

TEST_F(VerletClusterTowerTest, testCollectAllActualParticles) {
  ClusterTower<Particle> tower(4);
  constexpr size_t numParticles = 21;

  for (size_t i = 0; i < numParticles; i++) {
    tower.addParticle({{}, {}, i + 1});
  }
  tower.generateClusters();
  tower.setDummyValues(0, 0);

  const auto &actualParticles = tower.collectAllActualParticles();
  EXPECT_EQ(actualParticles.size(), numParticles);

  for (size_t i = 0; i < numParticles; i++) {
    // All particles have an ID not zero, only dummies have ID 0.
    EXPECT_NE(actualParticles[i].getID(), 0);
  }
}

TEST_F(VerletClusterTowerTest, testIterator) {
  constexpr size_t clusterSize = 4;
  ClusterTower<Particle> tower(clusterSize);
  constexpr size_t numParticles = 21;
  constexpr size_t numDummies = clusterSize - (numParticles % clusterSize);
  static_assert(numDummies != 0, "Without dummies this test is boring. Adapt parameters!");

  // map that counts how often each ID was found
  std::map<size_t, size_t> ids;
  for (size_t i = 0; i < numParticles; i++) {
    const auto id = i + 1;
    tower.addParticle({{}, {}, i + 1});
    ids.emplace(id, 0);
  }
  tower.generateClusters();
  tower.setDummyValues(0, 0);

  // Check that iterator iterates at least over particles, even dummies.
  size_t numDummiesFound = 0;
  for (const auto &particle : tower) {
    if (not(particle.isDummy())) {
      ++ids[particle.getID()];
      // if particle is NOT a dummy make sure their ID has not been encountered before
      EXPECT_EQ(ids[particle.getID()], 1)
          << "All particles should only be encountered exactly once! id: " << particle.getID();
    } else {
      ++numDummiesFound;
    }
  }
  // count found IDs in the map
  const size_t numParticlesFound =
      std::transform_reduce(ids.begin(), ids.end(), 0, std::plus(), [](const auto &pair) { return pair.second; });
  EXPECT_EQ(numParticlesFound, numParticles) << "Not all particles were found by the iterator.";
  EXPECT_EQ(numDummiesFound, numDummies) << "Not all dummies were found by the iterator.";
}