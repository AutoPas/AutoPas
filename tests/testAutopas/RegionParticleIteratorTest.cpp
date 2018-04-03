/**
 * @file RegionParticleIteratorTest.cpp
 * @author seckler
 * @date 03.04.18
 */

#include "RegionParticleIteratorTest.h"

using namespace autopas;

TEST_F(RegionParticleIteratorTest, testLinkedCellsRegionParticleIterator) {
  LinkedCells<TouchableParticle, FullParticleCell<TouchableParticle>>
      lcContainer(_boxMin, _boxMax, _cutoff);

  fillWithParticles(lcContainer);

  for (auto iterator = lcContainer.getRegionIterator(_regionMin, _regionMax);
       iterator.isValid(); ++iterator) {
    iterator->touch();
  }

  for (auto iterator = lcContainer.begin(); iterator.isValid(); ++iterator) {
    EXPECT_EQ(iterator->inBox(_regionMin, _regionMax) ? 1 : 0,
              iterator->getNumTouched());
//    std::cout << "id: " << iterator->getID() << " at [" << iterator->getR()[0]
//              << ", " << iterator->getR()[1] << ", " << iterator->getR()[2]
//              << "] touched:" << iterator->getNumTouched() << std::endl;
  }
}
