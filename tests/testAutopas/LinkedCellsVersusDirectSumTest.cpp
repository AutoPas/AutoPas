/*
 * LinkedCellsVersusDirectSumTest.cpp
 *
 *  Created on: 23 Jan 2018
 *      Author: tchipevn
 */

#include "LinkedCellsVersusDirectSumTest.h"
#include <cstdlib>

LinkedCellsVersusDirectSumTest::LinkedCellsVersusDirectSumTest()
    : _directSum(getBoxMin(), getBoxMax(), getCutoff()),
      _linkedCells(getBoxMin(), getBoxMax(), getCutoff()) {
  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  autopas::LJFunctor<autopas::MoleculeLJ>::setGlobals(getCutoff(), eps, sig,
                                                      shift);
}

double LinkedCellsVersusDirectSumTest::fRand(double fMin, double fMax) const {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

std::array<double, 3> LinkedCellsVersusDirectSumTest::randomPosition(
    const std::array<double, 3> &boxMin,
    const std::array<double, 3> &boxMax) const {
  std::array<double, 3> r{};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}

void LinkedCellsVersusDirectSumTest::fillContainerWithMolecules(
    unsigned long numMolecules,
    autopas::ParticleContainer<autopas::MoleculeLJ,
                               autopas::FullParticleCell<autopas::MoleculeLJ>>
        &cont) const {
  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(cont.getBoxMin()), boxMax(cont.getBoxMax());

  for (int i = 0; i < numMolecules; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::MoleculeLJ m(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    cont.addParticle(m);
  }
}

void LinkedCellsVersusDirectSumTest::test(unsigned long numMolecules,
                                          double rel_err_tolerance) {
  fillContainerWithMolecules(numMolecules, _directSum);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    _linkedCells.addParticle(*it);
  }

  autopas::LJFunctor<autopas::MoleculeLJ> func;
  _directSum.iteratePairwise2(&func);
  _linkedCells.iteratePairwise2(&func);

  auto itDirect = _directSum.begin();
  auto itLinked = _linkedCells.begin();

  std::vector<std::array<double, 3>> forcesDirect(numMolecules),
      forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
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

  autopas::FlopCounterFunctor<autopas::MoleculeLJ> flopsDirect(getCutoff()),
      flopsLinked(getCutoff());
  _directSum.iteratePairwise2(&flopsDirect);
  _linkedCells.iteratePairwise2(&flopsLinked);

  ASSERT_EQ(flopsLinked.getKernelCalls(), flopsDirect.getKernelCalls());
  ASSERT_LE(flopsLinked.getDistanceCalculations(),
            flopsDirect.getDistanceCalculations());
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
  double rel_err_tolerance = 1e-12;
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
