/**
 * @file LinkedCellsVersusDirectSumTest.cpp
 * @author tchipev
 * @date 23.01.18
 */

#include "LinkedCellsVersusDirectSumTest.h"
#include "testingHelpers/RandomGenerator.h"

LinkedCellsVersusDirectSumTest::LinkedCellsVersusDirectSumTest()
    : _directSum(getBoxMin(), getBoxMax(), getCutoff(), 0.),
      _linkedCells(getBoxMin(), getBoxMax(), getCutoff(), 0., 1.) {}

void LinkedCellsVersusDirectSumTest::test(unsigned long numMolecules, double rel_err_tolerance) {
  autopas::MoleculeLJ defaultParticle;
  RandomGenerator::fillWithParticles(_directSum, defaultParticle, numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    _linkedCells.addParticle(*it);
  }

  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  std::map<unsigned long, double> universalMap;
  for (unsigned long i = 0; i < numMolecules; i++) {
    universalMap.emplace(i, 1.0);
  }
  ParticleClassLibrary PCL = ParticleClassLibrary(universalMap, universalMap, universalMap);

  autopas::LJFunctor<Molecule, FMCell> func(getCutoff(), PCL, shift);

  autopas::C08Traversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos, true> traversalLJ(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &func);
  autopas::DirectSumTraversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos, true>
      traversalDS(&func);
  _directSum.iteratePairwise(&traversalDS);
  _linkedCells.iteratePairwise(&traversalLJ);

  auto itDirect = _directSum.begin();
  auto itLinked = _linkedCells.begin();

  std::vector<std::array<double, 3>> forcesDirect(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forcesDirect.at(m.getID()) = m.getF();
  }

  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forcesLinked.at(m.getID()) = m.getF();
  }

  for (size_t i = 0; i < numMolecules; ++i) {
    for (unsigned int d = 0; d < 3; ++d) {
      double f1 = forcesDirect[i][d];
      double f2 = forcesLinked[i][d];
      double abs_err = std::abs(f1 - f2);
      double rel_err = std::abs(abs_err / f1);
      EXPECT_LT(rel_err, rel_err_tolerance);
    }
  }

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsDirect(getCutoff()), flopsLinked(getCutoff());
  autopas::C08Traversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos, true>
      traversalFLOPSLC(_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &flopsLinked);
  autopas::DirectSumTraversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos,
                              true>
      traversalFLOPSDS(&flopsDirect);
  _directSum.iteratePairwise(&traversalFLOPSDS);
  _linkedCells.iteratePairwise(&traversalFLOPSLC);

  ASSERT_EQ(flopsLinked.getKernelCalls(), flopsDirect.getKernelCalls());
  ASSERT_LE(flopsLinked.getDistanceCalculations(), flopsDirect.getDistanceCalculations());
}

TEST_F(LinkedCellsVersusDirectSumTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-14;

  test(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusDirectSumTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusDirectSumTest, test2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;
  test(numMolecules, rel_err_tolerance);
}
