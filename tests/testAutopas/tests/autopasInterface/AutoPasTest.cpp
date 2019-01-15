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

/**
 * Check whether an AutoPas object can be rebuild using
 */
TEST_F(AutoPasTest, checkRebuilding) {

  autoPas.addParticle();
  // ensure some particles


  {
    decltype(autoPas) autoPasTmp;
    autoPasTmp.init({0., 0., 0.}, {10., 10., 10.}, 1., 0, 1, {autopas::ContainerOptions::linkedCells},
                    {autopas::TraversalOptions::c08});
    // ensure no particles


    // copy particles


    // move objects
    autoPas = std::move(autoPasTmp);
  }

  // ensure same particles as before
}