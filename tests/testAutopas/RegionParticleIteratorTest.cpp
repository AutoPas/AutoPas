/**
 * @file RegionParticleIteratorTest.cpp
 * @author seckler
 * @date 03.04.18
 */

#include "RegionParticleIteratorTest.h"

using namespace autopas;

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIterator) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0));

  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // check the touch. Iterating over cells provides more debug info than normal iterator.
  for (size_t cellId = 0; cellId < lcContainer.getCells().size(); ++cellId) {
    auto cellId3D = utils::ThreeDimensionalMapping::oneToThreeD(cellId, {7, 7, 7});
    for (auto pIter = lcContainer.getCells()[cellId].begin(); pIter.isValid(); ++pIter) {
      EXPECT_EQ(inBox(pIter->getR(), _regionMin, _regionMax) ? 1 : 0, pIter->getNumTouched())
          << "at: [" << pIter->getR()[0] << ", " << pIter->getR()[1] << ", " << pIter->getR()[2] << "]" << std::endl
          << "in cell: " << cellId << " [" << cellId3D[0] << " | " << cellId3D[1] << " | " << cellId3D[2] << "]";
    }
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorOwned) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0), 100);

  TouchableParticle part(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  lcContainer.addHaloParticle(part);

  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
                                                     autopas::IteratorBehavior::ownedOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    ASSERT_EQ(inBox(iterator->getR(), _boxMin, _regionMax) ? 1 : 0, iterator->getNumTouched());
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorBehaviorHalo) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0), 100);

  TouchableParticle part(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), 100);
  lcContainer.addHaloParticle(part);

  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax,
                                                     autopas::IteratorBehavior::haloOnly);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    ASSERT_EQ(inBox(iterator->getR(), ArrayMath::addScalar(_boxMin, -_cutoff * 0.5), _regionMax)
                  ? (inBox(iterator->getR(), _boxMin, _regionMax) ? 0 : 1)
                  : 0,
              iterator->getNumTouched());
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorEmpty) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add no particles

  int i = 0;
  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax); iterator.isValid(); ++iterator) {
    iterator->touch();
    i++;
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    ASSERT_EQ(inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched());
    i++;
  }
  ASSERT_EQ(i, 0);
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorCopyConstructor) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0));

  auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
  auto iterator2 = iterator;

  // touch them using the regionIterator
  for (; iterator2.isValid(); ++iterator2) {
    iterator2->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    ASSERT_EQ(inBox(iterator->getR(), _regionMin, _regionMax) ? 1 : 0, iterator->getNumTouched());
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorCopyAssignment) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer, TouchableParticle({0., 0., 0.}, 0));

  auto iterator2 = lcContainer.getRegionIterator(_regionMin, _regionMax);
  // touch them using the regionIterator
  for (; iterator2.isValid(); ++iterator2) {
    iterator2->touch();
  }

  auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
  // touch them using the regionIterator
  for (; iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  iterator2 = lcContainer.getRegionIterator(_regionMin, _regionMax);
  // touch them using the regionIterator
  for (; iterator2.isValid(); ++iterator2) {
    iterator2->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    ASSERT_EQ(inBox(iterator->getR(), _regionMin, _regionMax) ? 3 : 0, iterator->getNumTouched());
  }
}

/**
 * Tests for correct iterator behavior when some of the cells in and outside of the region of interest are empty
 */
TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorSparseDomain) {
  // box goes from {0,0,0} to {5,5,5} + one halo layer
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>> lcContainer(_boxMin, _boxMax, _cutoff);

  size_t idShouldTouch = 0;
  TouchableParticle p({0, 0, 0}, idShouldTouch);
  p.setR({2, 2, 2});
  p.setID(idShouldTouch++);
  lcContainer.addParticle(p);
  p.setR({2, 3, 2});
  p.setID(idShouldTouch++);
  lcContainer.addParticle(p);
  p.setR({2, 3, 3});
  p.setID(idShouldTouch++);
  lcContainer.addParticle(p);

  size_t idOffset = 1000;
  size_t idShouldNotTouch = idOffset;
  p.setR({1, 1, 1});
  p.setID(idShouldNotTouch++);
  lcContainer.addParticle(p);
  p.setR({2, 4.5, 2});
  p.setID(idShouldNotTouch++);
  lcContainer.addParticle(p);
  p.setR({4, 4, 4});
  p.setID(idShouldNotTouch++);
  lcContainer.addParticle(p);

  std::array<double, 3> regionOfInterstMin = {2, 2, 2};
  std::array<double, 3> regionOfInterstMax = {3, 4, 4};

  int particlesTouched = 0;
  for (auto iterator = lcContainer.getRegionIterator(regionOfInterstMin, regionOfInterstMax); iterator.isValid();
       ++iterator) {
    iterator->touch();
    ++particlesTouched;
  }
  EXPECT_EQ(particlesTouched, idShouldTouch);

  int particlesChecked = 0;
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    EXPECT_EQ(inBox(iterator->getR(), regionOfInterstMin, regionOfInterstMax) ? 1 : 0, iterator->getNumTouched());
    ++particlesChecked;
  }
  EXPECT_EQ(particlesChecked, idShouldTouch + idShouldNotTouch - idOffset);
}
