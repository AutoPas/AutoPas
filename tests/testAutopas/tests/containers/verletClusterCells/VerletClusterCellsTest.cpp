/**
 * @file VerletClusterCellsTest.cpp
 * @author jspahl
 * @date 6.4.19
 */

#include "VerletClusterCellsTest.h"

#include "autopas/containers/verletClusterCells/VerletClusterCells.h"
#include "autopas/containers/verletClusterCells/traversals/VerletClusterCellsTraversal.h"
#include "testingHelpers/TouchableParticle.h"

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

TEST_F(VerletClusterCellsTest, testVerletClusterBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 500);

  EXPECT_EQ(verletLists.updateContainer().size(), 0);
}

TEST_F(VerletClusterCellsTest, testNeighborListBuild) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 500);

  MockFunctor<Particle, FPCell> emptyFunctor;
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> dummyTraversal(
      &emptyFunctor, verletLists.getTraversalSelectorInfo().clusterSize);
  verletLists.rebuildNeighborLists(&dummyTraversal);
}

TEST_F(VerletClusterCellsTest, testVerletListIterator) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 64;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 500);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);
  std::vector<int> particlesOwn(500, 0);
  std::vector<int> particlesHalo(50, 0);
  std::vector<int> particlesBoth(500, 0);

  MockFunctor<Particle, FPCell> emptyFunctor;
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false>
      verletClusterCellsTraversal(&emptyFunctor, verletLists.getTraversalSelectorInfo().clusterSize);
  verletLists.rebuildNeighborLists(&verletClusterCellsTraversal);

  int numOwn = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    EXPECT_TRUE(iter->getID() < 500);
    ++particlesOwn[iter->getID()];
    ++numOwn;
  }
  EXPECT_GT(numOwn, 499);
  EXPECT_EQ(numOwn, 500);

  int numHalo = 0;

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    EXPECT_TRUE(iter->getID() < 50);
    ++particlesHalo[iter->getID()];
    ++numHalo;
  }
  EXPECT_GT(numHalo, 49);
  EXPECT_EQ(numHalo, 50);

  int numBoth = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    EXPECT_TRUE(iter->getID() < 500);
    ++particlesBoth[iter->getID()];
    ++numBoth;
  }
  EXPECT_EQ(numBoth, numOwn + numHalo);
  EXPECT_EQ(verletLists.getNumParticles(), 550);

  for (auto &it : particlesOwn) {
    EXPECT_EQ(it, 1);
  }
  for (auto &it : particlesHalo) {
    EXPECT_EQ(it, 1);
  }
  int i = 0;
  for (; i < 50; ++i) {
    EXPECT_EQ(particlesBoth[i], 2) << "on index " << i << std::endl;
  }
  for (; i < 500; ++i) {
    EXPECT_EQ(particlesBoth[i], 1) << "on index " << i << std::endl;
  }

  verletLists.deleteHaloParticles();

  EXPECT_FALSE(verletLists.begin(autopas::IteratorBehavior::haloOnly).isValid());
  // trigger rebuild
  verletLists.rebuildNeighborLists(&verletClusterCellsTraversal);
  std::vector<std::array<double, 3>> pos(500);

  i = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    EXPECT_TRUE(iter->getID() < 500);
    ++particlesBoth[iter->getID()];
    ++i;
  }
  EXPECT_EQ(i, 500);
  i = 0;
  for (; i < 50; ++i) {
    EXPECT_EQ(particlesBoth[i], 3) << "on index " << i << std::endl;
  }
  for (; i < 500; ++i) {
    EXPECT_EQ(particlesBoth[i], 2) << "on index " << i << std::endl;
  }
}

TEST_F(VerletClusterCellsTest, testVerletListIteratorDelete) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 64;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 500);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);

  std::vector<int> particlesBoth(500, 0);

  MockFunctor<Particle, FPCell> emptyFunctor;
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false>
      verletClusterCellsTraversal(&emptyFunctor, verletLists.getTraversalSelectorInfo().clusterSize);
  verletLists.rebuildNeighborLists(&verletClusterCellsTraversal);

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    if (iter->getID() % 2 == 0) {
      autopas::internal::deleteParticle(iter);
    }
  }

  int numBoth = 0;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    ++particlesBoth[iter->getID()];
    ++numBoth;
  }
  EXPECT_EQ(numBoth, 550 / 2);
  EXPECT_EQ(verletLists.getNumParticles(), 550 / 2);

  int i = 0;
  for (; i < 50; ++i) {
    EXPECT_EQ(particlesBoth[i], ((i % 2) == 0) ? 0 : 2) << "on index " << i << std::endl;
  }
  for (; i < 500; ++i) {
    EXPECT_EQ(particlesBoth[i], ((i % 2) == 0) ? 0 : 1) << "on index " << i << std::endl;
  }
}

TEST_F(VerletClusterCellsTest, testVerletParticleLoss) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 32;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 500);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);
  std::vector<int> particlesOwn(500, 0);
  std::vector<int> particlesHalo(50, 0);
  std::vector<int> particlesBoth(500, 0);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, false)).Times(AtLeast(1));
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false>
      verletClusterCellsTraversal(&emptyFunctor, verletLists.getTraversalSelectorInfo().clusterSize);

  verletLists.rebuildNeighborLists(&verletClusterCellsTraversal);

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

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    iter->setR(autopasTools::generators::RandomGenerator::randomPosition(min, max));
  }
  verletLists.iteratePairwise(&verletClusterCellsTraversal);
  verletLists.iteratePairwise(&verletClusterCellsTraversal);
  verletLists.iteratePairwise(&verletClusterCellsTraversal);

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

  for (auto &it : particlesOwn) {
    EXPECT_EQ(it, 2);
  }
  for (auto &it : particlesHalo) {
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

TEST_F(VerletClusterCellsTest, testIteratePairwiseWithoutNeighborlistRebuildThrows) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 32;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 50);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);

  MockFunctor<Particle, FPCell> emptyFunctor;
  EXPECT_CALL(emptyFunctor, AoSFunctor(_, _, false)).Times(0);
  autopas::VerletClusterCellsTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false>
      verletClusterCellsTraversal(&emptyFunctor, verletLists.getTraversalSelectorInfo().clusterSize);

  EXPECT_ANY_THROW(verletLists.iteratePairwise(&verletClusterCellsTraversal));
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

TEST_F(VerletClusterCellsTest, testUpdateHaloParticle) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {3, 3, 3};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 64;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 500);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);

  std::vector<Particle> haloToUpdate;
  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    if (iter->getID() % 3 == 0) {
      haloToUpdate.push_back(*iter);
    }
  }
  for (auto &it : haloToUpdate) {
    it.setF({1, 2, 3});
  }
  for (auto &it : haloToUpdate) {
    EXPECT_TRUE(verletLists.updateHaloParticle(it));
  }

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    if (iter->getID() % 3 == 0) {
      EXPECT_EQ(iter->getF()[0], 1);
      EXPECT_EQ(iter->getF()[1], 2);
      EXPECT_EQ(iter->getF()[2], 3);
    } else {
      EXPECT_EQ(iter->getF()[0], 0);
      EXPECT_EQ(iter->getF()[1], 0);
      EXPECT_EQ(iter->getF()[2], 0);
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

TEST_F(VerletClusterCellsTest, testVerletListRegionIterator) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {8, 8, 8};
  double cutoff = 1.;
  double skin = 0.2;
  int clusterSize = 64;
  autopas::VerletClusterCells<Particle> verletLists(min, max, cutoff, skin, clusterSize);

  autopasTools::generators::RandomGenerator::fillWithParticles(verletLists, Particle(), verletLists.getBoxMin(),
                                                               verletLists.getBoxMax(), 500);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(verletLists, Particle(), cutoff, 50);

  std::vector<int> particlesOwn(500, 0);
  std::vector<int> particlesHalo(50, 0);

  std::array<double, 3> minRegion = {3, 3, 3};
  std::array<double, 3> maxRegion = {7, 7, 5};

  for (auto iter = verletLists.getRegionIterator(minRegion, maxRegion, autopas::IteratorBehavior::haloAndOwned);
       iter.isValid(); ++iter) {
    EXPECT_TRUE(iter->getID() < 500);
    ++particlesOwn[iter->getID()];
  }

  EXPECT_FALSE(verletLists.getRegionIterator(minRegion, maxRegion, autopas::IteratorBehavior::haloOnly).isValid());

  for (auto iter = verletLists.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    if (particlesOwn[iter->getID()] > 0)
      EXPECT_TRUE(autopas::utils::inBox(iter->getR(), minRegion, maxRegion))
          << "On ID: " << iter->getID() << " position: (" << iter->getR()[0] << ", " << iter->getR()[1] << ", "
          << iter->getR()[2] << ")" << std::endl;
    else
      EXPECT_FALSE(autopas::utils::inBox(iter->getR(), minRegion, maxRegion))
          << "On ID: " << iter->getID() << " position: (" << iter->getR()[0] << ", " << iter->getR()[1] << ", "
          << iter->getR()[2] << ")" << std::endl;
  }
}