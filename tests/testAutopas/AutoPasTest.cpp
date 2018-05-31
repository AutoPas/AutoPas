/**
 * @file AutoPasTest.cpp
 * @author seckler
 * @date 29.05.18
 */

#include "AutoPasTest.h"
#include "AutoPas.h"

TEST_F(AutoPasTest, getRegionParticleIterator) {
  auto iter = autoPas.getRegionIterator({0., 0., 0.}, {4., 4., 4.});

  for (; iter.isValid(); ++iter) {
    iter->setR(iter->getR());
  }
}