/**
 * @file VerletClusterCellsTest.cpp
 * @author jspahl
 * @date 6.4.19
 */

#include "VerletClusterCellsTest.h"
#include "autopas/containers/verletClusterLists/VerletClusterCells.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClusterCellsTraversal.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

TEST_F(VerletClusterCellsTest, VerletListConstructor) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin);
}

TEST_F(VerletClusterCellsTest, testVerletListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletLists.addParticle(p);
  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletLists.addParticle(p2);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, false)).Times(AtLeast(1));
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> dummyTraversal(
      &emptyFunctor);
  verletLists.iteratePairwise(&emptyFunctor, &dummyTraversal, false);
}

TEST_F(VerletClusterCellsTest, testVerletListIterator) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 64;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, 2, clusterSize);

  RandomGenerator::fillWithParticles(verletLists, Particle(), 500);
  RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);
  std::vector<int> particlesOwn(500, 0);
  std::vector<int> particlesHalo(50, 0);
  std::vector<int> particlesBoth(500, 0);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, false)).Times(AtLeast(1));
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false>
      verletClusterCellsTraversal(&emptyFunctor);
  verletLists.iteratePairwise(&emptyFunctor, &verletClusterCellsTraversal, false);

  int numOwn = 0;
  int numDummyOwn = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    if (iter->getID() < 500)
      ++particlesOwn[iter->getID()];
    else
      ++numDummyOwn;
    ++numOwn;
  }
  EXPECT_GT(numOwn, 499);
  EXPECT_EQ(numOwn - numDummyOwn, 500);
  EXPECT_TRUE(numOwn % clusterSize == 0);

  int numHalo = 0;
  int numDummyHalo = 0;

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    if (iter->getID() < 50)
      ++particlesHalo[iter->getID()];
    else
      ++numDummyHalo;
    ++numHalo;
  }
  EXPECT_GT(numHalo, 49);
  EXPECT_EQ(numHalo - numDummyHalo, 50);
  EXPECT_TRUE(numHalo % clusterSize == 0);

  int numBoth = 0;
  int numDummyBoth = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    if (iter->getID() < 500)
      ++particlesBoth[iter->getID()];
    else
      ++numDummyBoth;
    ++numBoth;
  }
  EXPECT_EQ(numBoth - numDummyBoth, 550);
  EXPECT_EQ(numBoth, numOwn + numHalo);

  for (auto& it : particlesOwn) {
    EXPECT_EQ(it, 1);
  }
  for (auto& it : particlesHalo) {
    EXPECT_EQ(it, 1);
  }
  int i = 0;
  for (; i < 50; ++i) {
    EXPECT_EQ(particlesBoth[i], 2) << "on index " << i << std::endl;
  }
  for (; i < 500; ++i) {
    EXPECT_EQ(particlesBoth[i], 1) << "on index " << i << std::endl;
  }
}

TEST_F(VerletClusterCellsTest, testVerletParticleLoss) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 32;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, 2, clusterSize);

  RandomGenerator::fillWithParticles(verletLists, Particle(), 500);
  RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);
  std::vector<int> particlesOwn(500, 0);
  std::vector<int> particlesHalo(50, 0);
  std::vector<int> particlesBoth(500, 0);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, false)).Times(AtLeast(1));
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false>
      verletClusterCellsTraversal(&emptyFunctor);
  verletLists.iteratePairwise(&emptyFunctor, &verletClusterCellsTraversal, false);

  int numOwn = 0;
  int numDummyOwn = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    if (iter->getID() < 500)
      ++particlesOwn[iter->getID()];
    else
      ++numDummyOwn;

    ++numOwn;
  }
  EXPECT_GT(numOwn, 499);
  EXPECT_EQ(numOwn - numDummyOwn, 500);
  EXPECT_TRUE(numOwn % clusterSize == 0);

  int numHalo = 0;
  int numDummyHalo = 0;

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    if (iter->getID() < 50)
      ++particlesHalo[iter->getID()];
    else
      ++numDummyHalo;
    ++numHalo;
  }
  EXPECT_GT(numHalo, 49);
  EXPECT_EQ(numHalo - numDummyHalo, 50);
  EXPECT_TRUE(numHalo % clusterSize == 0);

  int numBoth = 0;
  int numDummyBoth = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    if (iter->getID() < 500)

      ++particlesBoth[iter->getID()];
    else
      ++numDummyBoth;
    ++numBoth;
  }
  EXPECT_EQ(numBoth - numDummyBoth, 550);
  EXPECT_EQ(numBoth, numOwn + numHalo);

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    iter->setR(RandomGenerator::randomPosition(min, max));
  }
  verletLists.iteratePairwise(&emptyFunctor, &verletClusterCellsTraversal, false);
  verletLists.iteratePairwise(&emptyFunctor, &verletClusterCellsTraversal, false);
  verletLists.iteratePairwise(&emptyFunctor, &verletClusterCellsTraversal, false);

  numOwn = 0;
  numDummyOwn = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    if (iter->getID() < 500)
      ++particlesOwn[iter->getID()];
    else
      ++numDummyOwn;

    ++numOwn;
  }
  EXPECT_GT(numOwn, 499);
  EXPECT_EQ(numOwn - numDummyOwn, 500);
  EXPECT_TRUE(numOwn % clusterSize == 0);

  numHalo = 0;
  numDummyHalo = 0;

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    if (iter->getID() < 50)
      ++particlesHalo[iter->getID()];
    else
      ++numDummyHalo;
    ++numHalo;
  }
  EXPECT_GT(numHalo, 49);
  EXPECT_EQ(numHalo - numDummyHalo, 50);
  EXPECT_TRUE(numHalo % clusterSize == 0);

  numBoth = 0;
  numDummyBoth = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    if (iter->getID() < 500)

      ++particlesBoth[iter->getID()];
    else
      ++numDummyBoth;
    ++numBoth;
  }
  EXPECT_EQ(numBoth - numDummyBoth, 550);
  EXPECT_EQ(numBoth, numOwn + numHalo);

  for (auto& it : particlesOwn) {
    EXPECT_EQ(it, 2);
  }
  for (auto& it : particlesHalo) {
    EXPECT_EQ(it, 2);
  }
  int i = 0;
  for (; i < 50; ++i) {
    EXPECT_EQ(particlesBoth[i], 4) << "on index " << i << std::endl;
  }
  for (; i < 500; ++i) {
    EXPECT_EQ(particlesBoth[i], 2) << "on index " << i << std::endl;
  }
}
