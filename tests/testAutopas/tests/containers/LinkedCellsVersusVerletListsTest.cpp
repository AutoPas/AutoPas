/**
 * @file LinkedCellsVersusDirectSumTest.cpp
 * @author seckler
 * @date 21.05.18
 */

#include "LinkedCellsVersusVerletListsTest.h"

LinkedCellsVersusVerletListsTest::LinkedCellsVersusVerletListsTest() : _verletLists(nullptr), _linkedCells(nullptr) {}

void LinkedCellsVersusVerletListsTest::test(unsigned long numMolecules, double rel_err_tolerance,
                                            std::array<double, 3> boxMax, bool useSoA, bool blackBoxMode) {
  _linkedCells = std::make_unique<lctype>(getBoxMin(), boxMax, getCutoff());
  if (blackBoxMode) {
    _verletLists = std::make_unique<vltype>(getBoxMin(), boxMax, getCutoff(), 0.1 * getCutoff(), 2,
                                            vltype::BuildVerletListType::VerletSoA, true);
  } else {
    _verletLists = std::make_unique<vltype>(getBoxMin(), boxMax, getCutoff(), 0.1 * getCutoff(), 2);
  }
  RandomGenerator::fillWithParticles(*_verletLists, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0), numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _verletLists->begin(); it.isValid(); ++it) {
    _linkedCells->addParticle(*it);
  }

  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  autopas::LJFunctor<Molecule, FMCell> func(getCutoff(), eps, sig, shift);

  autopas::C08Traversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, false, true> traversalLJ(
      _linkedCells->getCellBlock().getCellsPerDimensionWithHalo(), &func);
  _linkedCells->iteratePairwiseAoS(&func, &traversalLJ);

  autopas::DummyTraversal<FMCell> dt({0, 0, 0});  // The dummy doesn't need correct size of cellblock!
  if (useSoA) {
    _verletLists->iteratePairwiseSoA(&func, &dt);
  } else {
    _verletLists->iteratePairwiseAoS(&func, &dt);
  }

  auto itDirect = _verletLists->begin();
  auto itLinked = _linkedCells->begin();

  std::vector<std::array<double, 3>> forcesVerlet(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _verletLists->begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forcesVerlet.at(m.getID()) = m.getF();
  }

  for (auto it = _linkedCells->begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forcesLinked.at(m.getID()) = m.getF();
  }

  for (unsigned long i = 0; i < numMolecules; ++i) {
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesVerlet[i][d];
      double f2 = forcesLinked[i][d];
      EXPECT_NEAR(f1, f2, std::fabs(f1 * rel_err_tolerance));
    }
  }

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsVerlet(getCutoff()), flopsLinked(getCutoff());
  autopas::C08Traversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, false, true> traversalFLOPS(
      _linkedCells->getCellBlock().getCellsPerDimensionWithHalo(), &flopsLinked);
  _verletLists->iteratePairwiseAoS(&flopsVerlet, &traversalFLOPS);
  _linkedCells->iteratePairwiseAoS(&flopsLinked, &traversalFLOPS);

  EXPECT_EQ(flopsLinked.getKernelCalls(), flopsVerlet.getKernelCalls());

  // blackbox mode: the following line is only true, if the verlet lists do NOT use less cells than the linked cells
  // (for small scenarios), as the verlet lists fall back to linked cells.
  EXPECT_GE(flopsLinked.getDistanceCalculations(), flopsVerlet.getDistanceCalculations());
}

TEST_F(LinkedCellsVersusVerletListsTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-14;
  for (auto blackBox : {true, false}) {
    for (auto useSoA : {true, false}) {
      for (auto boxMax : {std::array<double, 3>{3., 3., 3.}, std::array<double, 3>{10., 10., 10.}}) {
        test(numMolecules, rel_err_tolerance, boxMax, useSoA, blackBox);
      }
    }
  }
}

TEST_F(LinkedCellsVersusVerletListsTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 2e-12;
  for (auto blackBox : {true, false}) {
    for (auto useSoA : {true, false}) {
      for (auto boxMax : {std::array<double, 3>{3., 3., 3.}, std::array<double, 3>{10., 10., 10.}}) {
        test(numMolecules, rel_err_tolerance, boxMax, useSoA, blackBox);
      }
    }
  }
}

TEST_F(LinkedCellsVersusVerletListsTest, test2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;
  for (auto blackBox : {true, false}) {
    for (auto useSoA : {true, false}) {
      for (auto boxMax : {std::array<double, 3>{3., 3., 3.}, std::array<double, 3>{10., 10., 10.}}) {
        test(numMolecules, rel_err_tolerance, boxMax, useSoA, blackBox);
      }
    }
  }
}