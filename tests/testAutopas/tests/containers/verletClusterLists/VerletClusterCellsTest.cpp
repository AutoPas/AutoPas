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

TEST_F(VerletClusterCellsTest, testParticleAdding) {
  autopas::VerletClusterCells<Particle> verletClusterCells({0., 0., 0.}, {10., 10., 10.}, 1.);
  int id = 1;
  for (double x : {-.5, 0., 5., 9.999, 10., 10.5}) {
    for (double y : {-.5, 0., 5., 9.999, 10., 10.5}) {
      for (double z : {-.5, 0., 5., 9.999, 10., 10.5}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
            z == 10.5) {
          EXPECT_ANY_THROW(verletClusterCells.addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(verletClusterCells.addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(verletClusterCells.addParticle(p));       // inside, therefore ok!
          EXPECT_ANY_THROW(verletClusterCells.addHaloParticle(p));  // inside, therefore not ok!
        }
      }
    }
  }
}

TEST_F(VerletClusterCellsTest, testGetNumParticles) {
  autopas::VerletClusterCells<Particle> verletClusterCells({0., 0., 0.}, {10., 10., 10.}, 1.);
  EXPECT_EQ(verletClusterCells.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletClusterCells.addParticle(p);
  EXPECT_EQ(verletClusterCells.getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletClusterCells.addParticle(p2);
  EXPECT_EQ(verletClusterCells.getNumParticles(), 2);
}

TEST_F(VerletClusterCellsTest, testDeleteAllParticles) {
  autopas::VerletClusterCells<Particle> verletClusterCells({0., 0., 0.}, {10., 10., 10.}, 1.);
  EXPECT_EQ(verletClusterCells.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  verletClusterCells.addParticle(p);
  EXPECT_EQ(verletClusterCells.getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  verletClusterCells.addParticle(p2);
  EXPECT_EQ(verletClusterCells.getNumParticles(), 2);

  verletClusterCells.deleteAllParticles();
  EXPECT_EQ(verletClusterCells.getNumParticles(), 0);
}

TEST_F(VerletClusterCellsTest, testCheckUpdateContainerNeededNoMove) {
  {
    autopas::VerletClusterCells<Particle> verletClusterCells({0., 0., 0.}, {10., 10., 10.}, 1.);
    int id = 1;
    for (double x : {-.5, 0., 5., 9.999, 10., 10.5}) {
      for (double y : {-.5, 0., 5., 9.999, 10., 10.5}) {
        for (double z : {-.5, 0., 5., 9.999, 10., 10.5}) {
          autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
          bool halo = false;
          for (int d = 0; d < 3; d++) {
            if (p.getR()[d] < 0. or p.getR()[d] >= 10.) {
              halo = true;
            }
          }
          if (halo) {
            verletClusterCells.addHaloParticle(p);
          } else {
            verletClusterCells.addParticle(p);
          }
          EXPECT_FALSE(verletClusterCells.isContainerUpdateNeeded());
        }
      }
    }
  }
  {
    autopas::VerletClusterCells<Particle> verletClusterCells({0., 0., 0.}, {10., 10., 10.}, 3.);
    int id = 1;
    for (double x : {-1.5, -.5, 0., 1. / 3, 2. / 3, 10., 10.5, 11.5}) {
      for (double y : {-1.5, -.5, 0., 1. / 3, 2. / 3, 10., 10.5, 11.5}) {
        for (double z : {-1.5, -.5, 0., 1. / 3, 2. / 3, 10., 10.5, 11.5}) {
          autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
          bool halo = false;
          for (int d = 0; d < 3; d++) {
            if (p.getR()[d] < 0. or p.getR()[d] >= 10.) {
              halo = true;
            }
          }
          if (halo) {
            verletClusterCells.addHaloParticle(p);
          } else {
            verletClusterCells.addParticle(p);
          }
          EXPECT_FALSE(verletClusterCells.isContainerUpdateNeeded());
        }
      }
    }
  }
}

TEST_F(VerletClusterCellsTest, testIsContainerUpdateNeeded) {
  std::array<double, 3> boxMin{0, 0, 0};
  std::array<double, 3> boxMax{10, 10, 10};
  double cutoff = 1.;
  autopas::VerletClusterCells<Particle> container(boxMin, boxMax, cutoff, 0, 2, 32);

  EXPECT_FALSE(container.isContainerUpdateNeeded());

  Particle p({1, 1, 1}, {0, 0, 0}, 0);
  container.addParticle(p);
  EXPECT_FALSE(container.isContainerUpdateNeeded());

  // Particle moves to halo cell -> needs update
  container.begin()->setR({-1, -1, -1});
  EXPECT_TRUE(container.isContainerUpdateNeeded());

  // Particle still in halo cell, but no container update called -> still needs update
  container.begin()->setR({-1.1, -0.9, -1});
  EXPECT_TRUE(container.isContainerUpdateNeeded());
}
