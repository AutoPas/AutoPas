/**
 * @file LinkedCellsTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "LinkedCellsTest.h"

TYPED_TEST_SUITE_P(LinkedCellsTest);

TYPED_TEST_P(LinkedCellsTest, testUpdateContainer) {
  using LinkedCellsType = TypeParam;
  LinkedCellsType linkedCells({0., 0., 0.}, {3., 3., 3.}, 1., 0., 1.);

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

TYPED_TEST_P(LinkedCellsTest, testUpdateContainerCloseToBoundary) {
  int id = 1;
  for (double x : {0., 5., 9.999}) {
    for (double y : {0., 5., 9.999}) {
      for (double z : {0., 5., 9.999}) {
        autopas::Particle p({x, y, z}, {0., 0., 0.}, id++);
        EXPECT_NO_THROW(this->_linkedCells.addParticle(p));  // inside, therefore ok!
      }
    }
  }
  std::set<unsigned long> movedIDs;
  // we move particles that are close to the boundary to outside of the container and remember the id's we moved
  for (auto iter = this->_linkedCells.begin(); iter.isValid(); ++iter) {
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
  auto invalidParticles = this->_linkedCells.updateContainer();

  // the particles should no longer be in the inner cells!
  for (auto iter = this->_linkedCells.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    EXPECT_EQ(movedIDs.count(iter->getID()), 0);
  }

  // the particles should now be inside of invalidParticles vector!
  EXPECT_EQ(movedIDs.size(), invalidParticles.size());
  for (auto &particle : invalidParticles) {
    EXPECT_EQ(movedIDs.count(particle.getID()), 1);
  }
}

REGISTER_TYPED_TEST_SUITE_P(LinkedCellsTest, testUpdateContainer, testUpdateContainerCloseToBoundary);

using MyTypes = ::testing::Types<autopas::LinkedCells<Particle>, autopas::LinkedCellsReferences<Particle>>;
INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LinkedCellsTest, MyTypes);
