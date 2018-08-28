/**
 * @file LinkedCellsVersusDirectSumTest.cpp
 * @author tchipev
 * @date 23.01.18
 */

#include "LinkedCellsVersusDirectSumTest.h"

LinkedCellsVersusDirectSumTest::LinkedCellsVersusDirectSumTest()
    : _directSum(getBoxMin(), getBoxMax(), getCutoff()), _linkedCells(getBoxMin(), getBoxMax(), getCutoff()) {
  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  autopas::LJFunctor<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>::setGlobals(getCutoff(), eps,
                                                                                                      sig, shift);
}

double LinkedCellsVersusDirectSumTest::fRand(double fMin, double fMax) const {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

std::array<double, 3> LinkedCellsVersusDirectSumTest::randomPosition(const std::array<double, 3> &boxMin,
                                                                     const std::array<double, 3> &boxMax) const {
  std::array<double, 3> r{};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}

void LinkedCellsVersusDirectSumTest::fillContainerWithMolecules(
    unsigned long numMolecules,
    autopas::ParticleContainer<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> &cont) const {
  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(cont.getBoxMin()), boxMax(cont.getBoxMax());

  for (unsigned long i = 0; i < numMolecules; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::MoleculeLJ m(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    cont.addParticle(m);
  }
}

void LinkedCellsVersusDirectSumTest::test(unsigned long numMolecules, double rel_err_tolerance) {
  fillContainerWithMolecules(numMolecules, _directSum);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    _linkedCells.addParticle(*it);
  }

  autopas::LJFunctor<Molecule, FMCell> func;
  autopas::C08Traversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, false, true> traversalLJ(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &func);
  _directSum.iteratePairwiseAoS(&func, &traversalLJ);
  _linkedCells.iteratePairwiseAoS(&func, &traversalLJ);

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
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesDirect[i][d];
      double f2 = forcesLinked[i][d];
      double abs_err = std::abs(f1 - f2);
      double rel_err = std::abs(abs_err / f1);
      EXPECT_LT(rel_err, rel_err_tolerance);
    }
  }

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsDirect(getCutoff()), flopsLinked(getCutoff());
  autopas::C08Traversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, false, true> traversalFLOPS(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &flopsLinked);
  _directSum.iteratePairwiseAoS(&flopsDirect, &traversalFLOPS);
  _linkedCells.iteratePairwiseAoS(&flopsLinked, &traversalFLOPS);

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
