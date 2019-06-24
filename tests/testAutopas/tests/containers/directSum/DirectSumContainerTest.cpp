/**
 * @file DirectSumContainerTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "DirectSumContainerTest.h"

TEST_F(DirectSumContainerTest, testParticleAdding) {
  autopas::DirectSum<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> directSum(
      {0., 0., 0.}, {10., 10., 10.}, 1., 0.);
  int id = 1;
  for (double x : {-.5, 0., 5., 9.999, 10., 10.5}) {
    for (double y : {-.5, 0., 5., 9.999, 10., 10.5}) {
      for (double z : {-.5, 0., 5., 9.999, 10., 10.5}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
            z == 10.5) {
          EXPECT_ANY_THROW(directSum.addParticle(p));     // outside, therefore not ok!
          EXPECT_NO_THROW(directSum.addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(directSum.addParticle(p));      // inside, therefore ok!
          EXPECT_NO_THROW(directSum.addHaloParticle(p));  // inside, but ok, as we have inprecise boundaries!
        }
      }
    }
  }
}

TEST_F(DirectSumContainerTest, testGetNumParticles) {
  autopas::DirectSum<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> directSum(
      {0., 0., 0.}, {10., 10., 10.}, 1., 0.);
  EXPECT_EQ(directSum.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  directSum.addParticle(p);
  EXPECT_EQ(directSum.getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  directSum.addParticle(p2);
  EXPECT_EQ(directSum.getNumParticles(), 2);
}

TEST_F(DirectSumContainerTest, testDeleteAllParticles) {
  autopas::DirectSum<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> directSum(
      {0., 0., 0.}, {10., 10., 10.}, 1., 0.);
  EXPECT_EQ(directSum.getNumParticles(), 0);

  std::array<double, 3> r = {2, 2, 2};
  Particle p(r, {0., 0., 0.}, 0);
  directSum.addParticle(p);
  EXPECT_EQ(directSum.getNumParticles(), 1);

  std::array<double, 3> r2 = {1.5, 2, 2};
  Particle p2(r2, {0., 0., 0.}, 1);
  directSum.addParticle(p2);
  EXPECT_EQ(directSum.getNumParticles(), 2);

  directSum.deleteAllParticles();
  EXPECT_EQ(directSum.getNumParticles(), 0);
}

TEST_F(DirectSumContainerTest, testIsContainerUpdateNeeded) {
  std::array<double, 3> boxMin{0, 0, 0};
  std::array<double, 3> boxMax{10, 10, 10};
  double cutoff = 1.;
  autopas::DirectSum<Particle, FPCell> container(boxMin, boxMax, cutoff, 0.);

  EXPECT_FALSE(container.isContainerUpdateNeeded());

  Particle p({1, 1, 1}, {0, 0, 0}, 0);
  container.addParticle(p);
  EXPECT_FALSE(container.isContainerUpdateNeeded());

  // Particle moves within cell -> needs no update
  container.begin()->setR({2.5, 1, 1});
  EXPECT_FALSE(container.isContainerUpdateNeeded());

  // Particle moves to different cell -> needs update
  container.begin()->setR({-1, -1, -1});
  EXPECT_TRUE(container.isContainerUpdateNeeded());
}

TEST_F(DirectSumContainerTest, testUpdateContainerCloseToBoundary) {
  autopas::DirectSum<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> directSum(
      {0., 0., 0.}, {10., 10., 10.}, 1., 0.);
  int id = 1;
  for (double x : {0., 5., 9.999}) {
    for (double y : {0., 5., 9.999}) {
      for (double z : {0., 5., 9.999}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        EXPECT_NO_THROW(directSum.addParticle(p));  // inside, therefore ok!
      }
    }
  }
  std::set<unsigned long> movedIDs;
  // we move particles that are close to the boundary to outside of the container and remember the id's we moved
  for (auto iter = directSum.begin(); iter.isValid(); ++iter) {
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
  auto invalidParticles = directSum.updateContainer();

  // the particles should no longer be in the inner cells!
  for (auto iter = directSum.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    EXPECT_EQ(movedIDs.count(iter->getID()), 0);
  }

  // the particles should now be inside of invalidParticles vector!
  EXPECT_EQ(movedIDs.size(), invalidParticles.size());
  for (auto &particle : invalidParticles) {
    EXPECT_EQ(movedIDs.count(particle.getID()), 1);
  }
}

TEST_F(DirectSumContainerTest, testUpdateContainerHalo) {
  autopas::DirectSum<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> directSum({0., 0., 0.},
                                                                                                {3., 3., 3.}, 1., 0.);

  autopas::Particle p({-0.5, -0.5, -0.5}, {0, 0, 0}, 42);
  directSum.addHaloParticle(p);

  // update container, will delete halo particles
  auto invalidParticles = directSum.updateContainer();
  // no particle should be returned
  EXPECT_EQ(invalidParticles.size(), 0);

  // no particle should remain
  auto iter = directSum.begin();
  EXPECT_FALSE(iter.isValid());
}