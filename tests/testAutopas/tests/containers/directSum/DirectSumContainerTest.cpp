/**
 * @file DirectSumContainerTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "DirectSumContainerTest.h"

TEST_F(DirectSumContainerTest, testUpdateContainerCloseToBoundary) {
  autopas::DirectSum<autopas::Particle> directSum({0., 0., 0.}, {10., 10., 10.}, 1., 0.);
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
  for (auto iter = directSum.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    EXPECT_EQ(movedIDs.count(iter->getID()), 0);
  }

  // the particles should now be inside of invalidParticles vector!
  EXPECT_EQ(movedIDs.size(), invalidParticles.size());
  for (auto &particle : invalidParticles) {
    EXPECT_EQ(movedIDs.count(particle.getID()), 1);
  }
}