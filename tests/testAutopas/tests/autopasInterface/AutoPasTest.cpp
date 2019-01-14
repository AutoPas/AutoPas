/**
 * @file AutoPasTest.cpp
 * @author seckler
 * @date 29.05.18
 */

#include "AutoPasTest.h"

/**
 * test whether the RegionIterator can be generated and something is returned.
 * This mainly makes certain, that the specific parts of the code can be compiled.
 */
TEST_F(AutoPasTest, getRegionParticleIterator) {
  auto iter = autoPas.getRegionIterator({0., 0., 0.}, {4., 4., 4.});

  for (; iter.isValid(); ++iter) {
    iter->setR(iter->getR());
  }
}