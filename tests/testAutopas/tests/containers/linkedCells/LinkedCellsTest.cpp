/**
 * @file LinkedCellsTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "LinkedCellsTest.h"

TYPED_TEST_SUITE_P(LinkedCellsTest);

TYPED_TEST_P(LinkedCellsTest, testUpdateContainer) {
  decltype(this->_linkedCells) linkedCells({0., 0., 0.}, {3., 3., 3.}, 1., 0., 1.);

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

  this->checkParticleIDsInCells(linkedCells, {{31ul, {0}}, {62ul, {1, 2}}, {63ul, {3}}, {93ul, {4}}}, true, __LINE__);

  // new locations for particles
  linkedCells.getCells()[31].begin()->setR({1.5, 0.5, 0.5});
  linkedCells.getCells()[62].begin()->setR({2.5, 1.5, 0.5});
  linkedCells.getCells()[63].begin()->setR({-0.5, -0.5, -0.5});
  linkedCells.getCells()[93].begin()->setR({1.6, 0.5, 0.5});

  auto invalidParticles = linkedCells.updateContainer(this->_keepListsValid);

  ASSERT_EQ(invalidParticles.size(), 1);
  EXPECT_EQ(invalidParticles[0].getID(), 3);

  if (this->_keepListsValid) {
    // if the lists are kept valid, particles are NOT moved between cells!
    this->checkParticleIDsInCells(linkedCells, {{31ul, {0}}, {62ul, {1, 2}}, {93ul, {4}}}, true, __LINE__);
  } else {
    // if the lists are not kept valid, particles should be moved between cells, so update the cells!
    this->checkParticleIDsInCells(linkedCells, {{32ul, {0, 4}}, {38ul, {1}}, {62ul, {2}}},
                                  false /*here, we do not know the order!*/, __LINE__);
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
  auto invalidParticles = this->_linkedCells.updateContainer(this->_keepListsValid);

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

// Workaround for storing two types.
// Currently, clang produces bugs if one tries to store this using a tuple or a pair.
// This problem is described in P0641R2 and occurs if an explicitly defaulted constructor cannot be instantiated.
// This is fixed in c++20.
template <typename first, typename second>
struct two_values {
  using first_t = first;
  using second_t = second;
};

using LC_true = two_values<autopas::LinkedCells<Particle>, std::true_type>;
using LC_false = two_values<autopas::LinkedCells<Particle>, std::false_type>;
using LCRef_true = two_values<autopas::LinkedCellsReferences<Particle>, std::true_type>;
using LCRef_false = two_values<autopas::LinkedCellsReferences<Particle>, std::false_type>;

using MyTypes = ::testing::Types<LC_true, LC_false, LCRef_true, LCRef_false>;

/// @todo c++20: replace with:
// using MyTypes = ::testing::Types<std::tuple<autopas::LinkedCells<Particle>, std::true_type>,
//                                  std::tuple<autopas::LinkedCells<Particle>, std::false_type>,
//                                  std::tuple<autopas::LinkedCellsReferences<Particle>, std::true_type>,
//                                  std::tuple<autopas::LinkedCellsReferences<Particle>, std::false_type> >;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LinkedCellsTest, MyTypes);
