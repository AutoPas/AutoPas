/**
 * @file LinkedCellsTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "LinkedCellsTest.h"
#include <gmock/gmock-generated-matchers.h>

LinkedCellsTest::LinkedCellsTest() : _linkedCells({0., 0., 0.}, {10., 10., 10.}, 1., 0., 1.) {}

TEST_F(LinkedCellsTest, testParticleAdding) {
  int id = 1;
  for (double x : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
    for (double y : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
      for (double z : {-1.5, -.5, 0., 5., 9.999, 10., 10.5, 11.5}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        if (x == -1.5 or y == -1.5 or z == -1.5 or x == 11.5 or y == 11.5 or z == 11.5) {
          EXPECT_ANY_THROW(_linkedCells.addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(_linkedCells.addHaloParticle(p));  // much outside, still ok because it is ignored!
        } else if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
                   z == 10.5) {
          EXPECT_ANY_THROW(_linkedCells.addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(_linkedCells.addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(_linkedCells.addParticle(p));  // inside, therefore ok!
          EXPECT_NO_THROW(
              _linkedCells.addHaloParticle(p));  // inside, but still ok, as halo particle can be inside of the domain!
        }
      }
    }
  }
}

TEST_F(LinkedCellsTest, testGetNumParticles) {
  EXPECT_EQ(_linkedCells.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  _linkedCells.addParticle(p);
  EXPECT_EQ(_linkedCells.getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  _linkedCells.addParticle(p2);
  EXPECT_EQ(_linkedCells.getNumParticles(), 2);
}

TEST_F(LinkedCellsTest, testDeleteAllParticles) {
  EXPECT_EQ(_linkedCells.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  _linkedCells.addParticle(p);
  EXPECT_EQ(_linkedCells.getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  _linkedCells.addParticle(p2);
  EXPECT_EQ(_linkedCells.getNumParticles(), 2);

  _linkedCells.deleteAllParticles();
  EXPECT_EQ(_linkedCells.getNumParticles(), 0);
}

TEST_F(LinkedCellsTest, testCheckUpdateContainerNeededNoMove) {
  {
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
            _linkedCells.addHaloParticle(p);
          } else {
            _linkedCells.addParticle(p);
          }
          EXPECT_FALSE(_linkedCells.isContainerUpdateNeeded());
        }
      }
    }
  }
  {
    decltype(_linkedCells) linkedCells2({0., 0., 0.}, {10., 10., 10.}, 3., 0., 1.);
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
            linkedCells2.addHaloParticle(p);
          } else {
            linkedCells2.addParticle(p);
          }
          EXPECT_FALSE(linkedCells2.isContainerUpdateNeeded());
        }
      }
    }
  }
}

TEST_F(LinkedCellsTest, testIsContainerUpdateNeeded) {
  EXPECT_FALSE(_linkedCells.isContainerUpdateNeeded());

  Particle p({1, 1, 1}, {0, 0, 0}, 0);
  _linkedCells.addParticle(p);
  EXPECT_FALSE(_linkedCells.isContainerUpdateNeeded());

  // Particle moves to different cell -> needs update
  _linkedCells.begin()->setR({2.5, 1, 1});
  EXPECT_TRUE(_linkedCells.isContainerUpdateNeeded());

  // Particle moves to halo cell -> needs update
  _linkedCells.begin()->setR({-1, -1, -1});
  EXPECT_TRUE(_linkedCells.isContainerUpdateNeeded());
}

TEST_F(LinkedCellsTest, testUpdateContainer) {
  autopas::LinkedCells<autopas::FullParticleCell<autopas::Particle>> linkedCells({0., 0., 0.}, {3., 3., 3.}, 1., 0.,
                                                                                 1.);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);

  linkedCells.addParticle(p1);
  linkedCells.addParticle(p2);
  linkedCells.addParticle(p3);
  linkedCells.addParticle(p4);
  linkedCells.addParticle(p5);

  // check particles are where we expect them to be (and nothing else)
  for (size_t i = 0; i < linkedCells.getCells().size(); ++i) {
    if (i == 31) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 1);
      EXPECT_EQ(linkedCells.getCells()[i].begin()->getID(), 0);
    } else if (i == 62) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 2);
      EXPECT_EQ(linkedCells.getCells()[i].begin()->getID(), 1);
      EXPECT_EQ((++(linkedCells.getCells()[i].begin()))->getID(), 2);
    } else if (i == 63) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 1);
      EXPECT_EQ(linkedCells.getCells()[i].begin()->getID(), 3);
    } else if (i == 93) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 1);
      EXPECT_EQ(linkedCells.getCells()[i].begin()->getID(), 4);
    } else {
      EXPECT_FALSE(linkedCells.getCells()[i].isNotEmpty());
    }
  }

  // new locations for particles
  linkedCells.getCells()[31].begin()->setR({1.5, 0.5, 0.5});
  linkedCells.getCells()[62].begin()->setR({2.5, 1.5, 0.5});
  linkedCells.getCells()[63].begin()->setR({-0.5, -0.5, -0.5});
  linkedCells.getCells()[93].begin()->setR({1.6, 0.5, 0.5});
  auto invalidParticles = linkedCells.updateContainer();

  ASSERT_EQ(invalidParticles.size(), 1);
  EXPECT_EQ(invalidParticles[0].getID(), 3);

  // verify particles are in correct new cells (and nowhere else)
  for (size_t i = 0; i < linkedCells.getCells().size(); ++i) {
    if (i == 0) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 0);
    } else if (i == 32) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 2);
      auto pIter = linkedCells.getCells()[i].begin();
      auto ids = {pIter->getID(), (++pIter)->getID()};
      EXPECT_THAT(ids, testing::UnorderedElementsAre(0, 4));
    } else if (i == 38) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 1);
      EXPECT_EQ(linkedCells.getCells()[i].begin()->getID(), 1);
    } else if (i == 62) {
      EXPECT_EQ(linkedCells.getCells()[i].numParticles(), 1);
      EXPECT_EQ(linkedCells.getCells()[i].begin()->getID(), 2);
    } else {
      EXPECT_FALSE(linkedCells.getCells()[i].isNotEmpty());
    }
  }
}

TEST_F(LinkedCellsTest, testUpdateContainerCloseToBoundary) {
  int id = 1;
  for (double x : {0., 5., 9.999}) {
    for (double y : {0., 5., 9.999}) {
      for (double z : {0., 5., 9.999}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        EXPECT_NO_THROW(_linkedCells.addParticle(p));  // inside, therefore ok!
      }
    }
  }
  std::set<unsigned long> movedIDs;
  // we move particles that are close to the boundary to outside of the container and remember the id's we moved
  for (auto iter = _linkedCells.begin(); iter.isValid(); ++iter) {
    for (unsigned short dim = 0; dim < 3; ++dim) {
      if (iter->getR()[dim] < 0.5) {
        auto r = iter->getR();
        // smallest double smaller than 0
        r[dim] = std::nexttoward(0., -1.);
        iter->setR(r);
        movedIDs.insert(iter->getID());
      }
      if (iter->getR()[dim] > 9.5) {
        auto r = iter->getR();
        r[dim] = 10.;
        iter->setR(r);
        movedIDs.insert(iter->getID());
      }
    }
  }

  // now update the container!
  auto invalidParticles = _linkedCells.updateContainer();

  // the particles should no longer be in the inner cells!
  for (auto iter = _linkedCells.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    EXPECT_EQ(movedIDs.count(iter->getID()), 0);
  }

  // the particles should now be inside of invalidParticles vector!
  EXPECT_EQ(movedIDs.size(), invalidParticles.size());
  for (auto &particle : invalidParticles) {
    EXPECT_EQ(movedIDs.count(particle.getID()), 1);
  }
}

TEST_F(LinkedCellsTest, testUpdateContainerHalo) {
  autopas::LinkedCells<autopas::FullParticleCell<autopas::Particle>> linkedCells({0., 0., 0.}, {3., 3., 3.}, 1., 0.,
                                                                                 1.);

  autopas::Particle p({-0.5, -0.5, -0.5}, {0, 0, 0}, 42);
  linkedCells.addHaloParticle(p);

  EXPECT_EQ(linkedCells.getCells()[0].numParticles(), 1);
  EXPECT_EQ(linkedCells.getCells()[0].begin()->getID(), 42);

  auto invalidParticles = linkedCells.updateContainer();

  // no particle should be returned
  EXPECT_EQ(invalidParticles.size(), 0);

  // no particle should remain
  auto iter = linkedCells.begin();
  EXPECT_FALSE(iter.isValid());
}