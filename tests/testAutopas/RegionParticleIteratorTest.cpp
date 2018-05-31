/**
 * @file RegionParticleIteratorTest.cpp
 * @author seckler
 * @date 03.04.18
 */

#include "RegionParticleIteratorTest.h"
#include "testingHelpers/RandomGenerator.h"

using namespace autopas;

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIterator) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>>
      lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer,
                                     TouchableParticle({0., 0., 0.}, 0));

  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    ASSERT_EQ(iterator->inBox(_regionMin, _regionMax) ? 1 : 0,
              iterator->getNumTouched());
  }
}

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIteratorEmpty) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>>
      lcContainer(_boxMin, _boxMax, _cutoff);

  // add no particles

  int i = 0;
  // touch them using the regionIterator
  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
       iterator.isValid(); ++iterator) {
    iterator->touch();
    i++;
  }

  // check the touch using the normal iterator
  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    //  std::cout << "id: " << iterator->getID() << " at [" <<
    //  iterator->getR()[0]
    //         << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
    //              << "] touched:" << iterator->getNumTouched() << std::endl;

    ASSERT_EQ(iterator->inBox(_regionMin, _regionMax) ? 1 : 0,
              iterator->getNumTouched());
    i++;
  }
  ASSERT_EQ(i, 0);
}

TEST_F(RegionParticleIteratorTest,
       testLinkedCellsRegionParticleIteratorCopyConstructor) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>>
      lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer,
                                     TouchableParticle({0., 0., 0.}, 0));

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

    ASSERT_EQ(iterator->inBox(_regionMin, _regionMax) ? 1 : 0,
              iterator->getNumTouched());
  }
}

TEST_F(RegionParticleIteratorTest,
       testLinkedCellsRegionParticleIteratorCopyAssignment) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>>
      lcContainer(_boxMin, _boxMax, _cutoff);

  // add a number of particles
  RandomGenerator::fillWithParticles(lcContainer,
                                     TouchableParticle({0., 0., 0.}, 0));

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

    ASSERT_EQ(iterator->inBox(_regionMin, _regionMax) ? 3 : 0,
              iterator->getNumTouched());
  }
}
