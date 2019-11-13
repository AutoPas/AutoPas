/**
 * @file DSCudaTraversalVersusDirectSumTest.cpp
 * @author jspahl
 * @date 11.03.19
 */

#include "DSCudaTraversalVersusDirectSumTest.h"

DSCudaTraversalVersusDirectSumTest::DSCudaTraversalVersusDirectSumTest()
    : _directSum(getBoxMin(), getBoxMax(), getCutoff(), 0.),
      _directSumCuda(getBoxMin(), getBoxMax(), getCutoff(), 0.) {}

double DSCudaTraversalVersusDirectSumTest::fRand(double fMin, double fMax) const {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

std::array<double, 3> DSCudaTraversalVersusDirectSumTest::randomPosition(const std::array<double, 3> &boxMin,
                                                                         const std::array<double, 3> &boxMax) const {
  std::array<double, 3> r{};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}

void DSCudaTraversalVersusDirectSumTest::fillContainerWithMolecules(
    unsigned long numMolecules,
    autopas::ParticleContainer<autopas::FullParticleCell<autopas::MoleculeLJ>> &cont) const {
  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(cont.getBoxMin()), boxMax(cont.getBoxMax());

  for (unsigned long i = 0; i < numMolecules; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::MoleculeLJ m(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    cont.addParticle(m);
  }
}

template <bool useNewton3, bool calculateGlobals>
void DSCudaTraversalVersusDirectSumTest::test(unsigned long numMolecules, double rel_err_tolerance) {
  fillContainerWithMolecules(numMolecules, _directSum);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    _directSumCuda.addParticle(*it);
  }

  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  autopas::LJFunctor<Molecule, FMCell, autopas::FunctorN3Modes::Both, calculateGlobals> funcDS(getCutoff(), eps, sig,
                                                                                               shift);
  autopas::LJFunctor<Molecule, FMCell, autopas::FunctorN3Modes::Both, calculateGlobals> funcDScuda(getCutoff(), eps,
                                                                                                   sig, shift);

  autopas::DirectSumTraversal<FMCell,
                              autopas::LJFunctor<Molecule, FMCell, autopas::FunctorN3Modes::Both, calculateGlobals>,
                              autopas::DataLayoutOption::cuda, useNewton3>
      traversalDSCuda(&funcDScuda, getCutoff());
  autopas::DirectSumTraversal<FMCell,
                              autopas::LJFunctor<Molecule, FMCell, autopas::FunctorN3Modes::Both, calculateGlobals>,
                              autopas::DataLayoutOption::aos, useNewton3>
      traversalDS(&funcDS, getCutoff());

  funcDS.initTraversal();
  _directSum.iteratePairwise(&traversalDS);
  funcDS.endTraversal(useNewton3);

  funcDScuda.initTraversal();
  _directSumCuda.iteratePairwise(&traversalDSCuda);
  funcDScuda.endTraversal(useNewton3);

  auto itDirect = _directSum.begin();
  auto itLinked = _directSumCuda.begin();

  std::vector<std::array<double, 3>> forcesDirect(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _directSum.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    forcesDirect.at(m.getID()) = m.getF();
  }

  for (auto it = _directSumCuda.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
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
    EXPECT_LT(std::abs((funcDS.getUpot() - funcDScuda.getUpot()) / funcDS.getUpot()), rel_err_tolerance);
    EXPECT_LT(std::abs((funcDS.getVirial() - funcDScuda.getVirial()) / funcDS.getVirial()), rel_err_tolerance);
  }
}
#if defined(AUTOPAS_CUDA)
TEST_F(DSCudaTraversalVersusDirectSumTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<false>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, test500) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<false>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<false>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, testN3100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, testN3500) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, testN31000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, test100Globals) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<false, true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, test500Globals) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<false, true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, test1000Globals) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<false, true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, testN3100Globals) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-13;

  test<true, true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, testN3500Globals) {
  unsigned long numMolecules = 500;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-12;

  test<true, true>(numMolecules, rel_err_tolerance);
}

TEST_F(DSCudaTraversalVersusDirectSumTest, testN31000Globals) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-12;
  test<true, true>(numMolecules, rel_err_tolerance);
}
#endif
