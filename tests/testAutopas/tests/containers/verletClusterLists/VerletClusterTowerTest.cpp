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
    EXPECT_EQ(tower.getNumDummyParticles(), clusterSize - numParticlesLastCluster);

    EXPECT_EQ(tower.getNumClusters(), (numParticles + tower.getNumDummyParticles()) / clusterSize);

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
    for (size_t i = 1; i <= tower.getNumDummyParticles(); i++) {
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

  for (size_t i = 0; i < numParticles; i++) {
    tower.addParticle({{}, {}, i + 1});
  }
  tower.generateClusters();
  tower.setDummyValues(0, 0);

  // Check that iterator iterates over all these particles and not over dummies. Dummies would have ID 0.
  std::vector<size_t> IDs;
  IDs.reserve(numParticles);
  size_t numDummiesFound = 0;
  for (const auto &particle : tower) {
    if (particle.getID(), std::numeric_limits<size_t>::max()) {
      ++numDummiesFound;
    } else {
      EXPECT_TRUE(std::find(IDs.begin(), IDs.end(), particle.getID()) == IDs.end());
      IDs.push_back(particle.getID());
    }
  }
}