/**
 * @file LinkedCellsVersusVerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "LinkedCellsVersusVerletClusterListsTest.h"

LinkedCellsVersusVerletClusterListsTest::LinkedCellsVersusVerletClusterListsTest()
    : _verletLists(getBoxMin(), getBoxMax(), getCutoff(), 0.1 * getCutoff(), 2),
      _linkedCells(getBoxMin(), getBoxMax(), getCutoff()) {
  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  autopas::LJFunctor<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>::setGlobals(getCutoff(), eps,
                                                                                                      sig, shift);
}

void LinkedCellsVersusVerletClusterListsTest::test(unsigned long numMolecules, double rel_err_tolerance) {
  RandomGenerator::fillWithParticles(_linkedCells, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0), numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    _verletLists.addParticle(*it);
  }

  // TODO tolerance test

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsVerlet(getCutoff()), flopsLinked(getCutoff());
  autopas::C08Traversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, false, false> traversalFLOPS(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &flopsLinked);
  _verletLists.iteratePairwiseAoS(&flopsVerlet, &traversalFLOPS, false);
  _linkedCells.iteratePairwiseAoS(&flopsLinked, &traversalFLOPS, false);

  ASSERT_EQ(flopsLinked.getKernelCalls(), flopsVerlet.getKernelCalls());
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-14;

  test(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 2e-12;
  test(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, test2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;
  test(numMolecules, rel_err_tolerance);
}
