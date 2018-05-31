/*
 * LinkedCellsVersusDirectSumTest.cpp
 *
 *  Created on: 23 Jan 2018
 *      Author: tchipevn
 */

#include "LinkedCellsVersusVerletListsTest.h"
#include <testingHelpers/RandomGenerator.h>
#include <cstdlib>

LinkedCellsVersusVerletListsTest::LinkedCellsVersusVerletListsTest()
    : _verletLists(getBoxMin(), getBoxMax(), getCutoff(), 0.1 * getCutoff(), 2),
      _linkedCells(getBoxMin(), getBoxMax(), getCutoff()) {
  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  autopas::LJFunctor<
      autopas::MoleculeLJ,
      autopas::FullParticleCell<autopas::MoleculeLJ>>::setGlobals(getCutoff(),
                                                                  eps, sig,
                                                                  shift);
}

void LinkedCellsVersusVerletListsTest::test(unsigned long numMolecules,
                                            double rel_err_tolerance) {
  RandomGenerator::fillWithParticles(
      _verletLists, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0),
      numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _verletLists.begin(); it.isValid(); ++it) {
    _linkedCells.addParticle(*it);
  }

  autopas::LJFunctor<autopas::MoleculeLJ,
                     autopas::FullParticleCell<autopas::MoleculeLJ>>
      func;
  _verletLists.iteratePairwiseAoS2(&func);
  _linkedCells.iteratePairwiseAoS2(&func);

  auto itDirect = _verletLists.begin();
  auto itLinked = _linkedCells.begin();

  std::vector<std::array<double, 3>> forcesDirect(numMolecules),
      forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _verletLists.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forcesDirect.at(m.getID()) = m.getF();
  }

  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forcesLinked.at(m.getID()) = m.getF();
  }

  for (int i = 0; i < numMolecules; ++i) {
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesDirect[i][d];
      double f2 = forcesLinked[i][d];
      double abs_err = std::abs(f1 - f2);
      double rel_err = std::abs(abs_err / f1);
      EXPECT_LT(rel_err, rel_err_tolerance);
    }
  }

  autopas::FlopCounterFunctor<autopas::MoleculeLJ,
                              autopas::FullParticleCell<autopas::MoleculeLJ>>
      flopsVerlet(getCutoff()), flopsLinked(getCutoff());
  _verletLists.iteratePairwiseAoS2(&flopsVerlet);
  _linkedCells.iteratePairwiseAoS2(&flopsLinked);

  ASSERT_EQ(flopsLinked.getKernelCalls(), flopsVerlet.getKernelCalls());
  ASSERT_GE(flopsLinked.getDistanceCalculations(),
            flopsVerlet.getDistanceCalculations());
}

TEST_F(LinkedCellsVersusVerletListsTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-14;

  test(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusVerletListsTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;
  test(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusVerletListsTest, test2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;
  test(numMolecules, rel_err_tolerance);
}
