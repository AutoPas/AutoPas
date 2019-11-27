/**
 * @file LCCudaTraversalVersusDirectSumTest.cpp
 * @author jspahl
 * @date 11.03.19
 */

#include "LCCudaTraversalVersusDirectSumTest.h"

#include "autopas/containers/linkedCells/traversals/C01CudaTraversal.h"
//#include "autopas/pairwiseFunctors/LJFunctor.h"

LCCudaTraversalVersusDirectSumTest::LCCudaTraversalVersusDirectSumTest()
    : _directSum(getBoxMin(), getBoxMax(), getCutoff(), 0.),
      _linkedCells(getBoxMin(), getBoxMax(), getCutoff(), 0., 1. /*cell size factor*/) {}

double LCCudaTraversalVersusDirectSumTest::fRand(double fMin, double fMax) const {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

std::array<double, 3> LCCudaTraversalVersusDirectSumTest::randomPosition(const std::array<double, 3> &boxMin,
                                                                         const std::array<double, 3> &boxMax) const {
  std::array<double, 3> r{};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}

void LCCudaTraversalVersusDirectSumTest::fillContainerWithMolecules(unsigned long numMolecules,
                                                                    autopas::ParticleContainer<FMCell> &cont) const {
  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(cont.getBoxMin()), boxMax(cont.getBoxMax());

  for (unsigned long i = 0; i < numMolecules; ++i) {
    auto id = static_cast<unsigned long>(i);
    Molecule m(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    cont.addParticle(m);
  }
}

template <bool useNewton3, bool calculateGlobals>
void LCCudaTraversalVersusDirectSumTest::test(unsigned long numMolecules, double rel_err_tolerance) {
  fillContainerWithMolecules(numMolecules, _directSum);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    _linkedCells.addParticle(*it);
  }

  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::LJFunctor<Molecule, FMCell, /*mixing*/ false, autopas::FunctorN3Modes::Both, calculateGlobals> funcDS(
      getCutoff(), shift);
  funcDS.setParticleProperties(eps * 24, sig * sig);
  autopas::LJFunctor<Molecule, FMCell, /*mixing*/ false, autopas::FunctorN3Modes::Both, calculateGlobals> funcLC(
      getCutoff(), shift);
  funcLC.setParticleProperties(eps * 24, sig * sig);

  autopas::DirectSumTraversal<FMCell, decltype(funcDS), autopas::DataLayoutOption::aos, useNewton3> traversalDS(
      &funcDS, getCutoff());
  autopas::C01CudaTraversal<FMCell, decltype(funcLC), autopas::DataLayoutOption::cuda, useNewton3> traversalLC(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &funcLC);

  funcDS.initTraversal();
  _directSum.iteratePairwise(&traversalDS);
  funcDS.endTraversal(useNewton3);

  funcLC.initTraversal();
  _linkedCells.iteratePairwise(&traversalLC);
  funcLC.endTraversal(useNewton3);

  auto itDirect = _directSum.begin();
  auto itLinked = _linkedCells.begin();

  std::vector<std::array<double, 3>> forcesDirect(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    Molecule &m = *it;
    forcesDirect.at(m.getID()) = m.getF();
  }

  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    Molecule &m = *it;
    forcesLinked.at(m.getID()) = m.getF();
  }

  for (size_t i = 0; i < numMolecules; ++i) {
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesDirect[i][d];
      double f2 = forcesLinked[i][d];
      double abs_err = std::abs(f1 - f2);
      double rel_err = std::abs(abs_err / f1);
      EXPECT_LT(rel_err, rel_err_tolerance) << " for ParticleID: " << i << " dim:" << d << " " << f1 << " vs " << f2;
    }
  }

  if (calculateGlobals) {
    EXPECT_LT(std::abs((funcDS.getUpot() - funcLC.getUpot()) / funcDS.getUpot()), rel_err_tolerance);
    EXPECT_LT(std::abs((funcDS.getVirial() - funcLC.getVirial()) / funcDS.getVirial()), rel_err_tolerance);
  }
}
#if defined(AUTOPAS_CUDA)
TEST_F(LCCudaTraversalVersusDirectSumTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<false>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, test500) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<false>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<false>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, testN3100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, testN3500) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, testN31000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, test100Globals) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<false, true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, test500Globals) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<false, true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, test1000Globals) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<false, true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, testN3100Globals) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<true, true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, testN3500Globals) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<true, true>(numMolecules, rel_err_tolerance);
}

TEST_F(LCCudaTraversalVersusDirectSumTest, testN31000Globals) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<true, true>(numMolecules, rel_err_tolerance);
}
#endif
