/**
 * @file LJFunctorFlopCounterTest.cpp
 * @author F. Gratl
 * @date 01.06.18
 */

#include "LJFunctorFlopCounterTest.h"

#include "autopas/AutoPasDecl.h"
#include "autopas/utils/WrapOpenMP.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;
/**
 * Generates a square of four particles, iterates over it with the LJFunctor and checks the values of getNumFLOPs() and getHitRate()
 * @param dataLayoutOption
 * @param newton3
 */
void LJFunctorFlopCounterTest::testFLOPCounter(autopas::DataLayoutOption dataLayoutOption, bool newton3) {
  autopas::AutoPas<Molecule> autoPas;

  autoPas.setBoxMin({0, 0, 0});
  autoPas.setBoxMax({3, 3, 3});
  autoPas.setCutoff(1);
  autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
  autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
  if (newton3) {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::enabled});
  } else {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  }

  autoPas.init();

  const std::vector<Molecule> molVec{Molecule({1, 1, 1}, {0, 0, 0}, 0), Molecule({1, 1, 2}, {0, 0, 0}, 1),
                               Molecule({1, 2, 1}, {0, 0, 0}, 2), Molecule({1, 2, 2}, {0, 0, 0}, 3)};

  for (auto &m : molVec) {
    autoPas.addParticle(m);
  }

  mdLib::LJFunctor<Molecule, false, false, autopas::FunctorN3Modes::Both, false, true, true> ljFunctor(autoPas.getCutoff());

  autoPas.iteratePairwise(&ljFunctor);

  // every particle checks the distance to all others. If newton3, only half of the calculations are made due to Newton 3.
  const auto expectedDistanceCalculations = newton3 ? molVec.size() * (molVec.size() - 1) / 2 : molVec.size() * (molVec.size() - 1);

  // in theory each particle has two in range but only one kernel call because of Newton 3.
  // Each particle has two others in range -> 2 kernel calls. With newton3, only half of these kernel calls happen.
  const auto expectedKernelCalls = newton3 ? molVec.size() : 2 * molVec.size();

  // distance calculations cost 8 flops, LJ kernel calls without Newton3 cost 15 FLOPs, with Newton 3 cost 18 flops
  const int numFLOPsPerKernelCall = newton3 ? 18 : 15;
  const auto expectedFlops = expectedDistanceCalculations * 8 + expectedKernelCalls * numFLOPsPerKernelCall;
  ASSERT_EQ(expectedFlops, ljFunctor.getNumFLOPs());

  // two out of three particles are in range
  const auto expectedHitRate = 2. / 3.;
  ASSERT_NEAR(expectedHitRate, ljFunctor.getHitRate(), 1e-14);
}

TEST_F(LJFunctorFlopCounterTest, testFlopCounterAoS4N3Mol) { testFLOPCounter(autopas::DataLayoutOption::aos, true); }

TEST_F(LJFunctorFlopCounterTest, testFlopCounterSoA4N3Mol) { testFLOPCounter(autopas::DataLayoutOption::soa, true); }

TEST_F(LJFunctorFlopCounterTest, testFlopCounterAoS4NoN3Mol) { testFLOPCounter(autopas::DataLayoutOption::aos, false); }

TEST_F(LJFunctorFlopCounterTest, testFlopCounterSoA4NoN3Mol) { testFLOPCounter(autopas::DataLayoutOption::soa, false); }

TEST_F(LJFunctorFlopCounterTest, testFlopCounterAoSN3OpenMP) {
  bool newton3 = true;
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 2., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 2.2, 0.3}, {0., 0., 0.}, 1, 0);

  double cutoff = 1.;

  mdLib::LJFunctor<Molecule> ljFunctor(cutoff);
  autopas::FlopCounterFunctor<Molecule, mdLib::LJFunctor<Molecule>> functor(ljFunctor, cutoff);

  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
  AUTOPAS_OPENMP(parallel) {
    AUTOPAS_OPENMP(sections) {
      AUTOPAS_OPENMP(section)
      functor.AoSFunctor(p1, p2, newton3);
      AUTOPAS_OPENMP(section)
      functor.AoSFunctor(p3, p4, newton3);
    }
  }
}

TEST_F(LJFunctorFlopCounterTest, testFlopCounterSoAOpenMP) {
  bool newton3 = true;
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p5({1., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p6({1.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p7({1., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p8({1.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  double cutoff = 1.;

  mdLib::LJFunctor<Molecule> ljFunctor(cutoff);
  autopas::FlopCounterFunctor<Molecule, mdLib::LJFunctor<Molecule>> functor(ljFunctor, cutoff);

  autopas::FullParticleCell<Molecule> cell1;
  cell1.addParticle(p1);
  cell1.addParticle(p2);

  autopas::FullParticleCell<Molecule> cell2;
  cell2.addParticle(p3);
  cell2.addParticle(p4);

  autopas::FullParticleCell<Molecule> cell3;
  cell3.addParticle(p5);
  cell3.addParticle(p6);

  autopas::FullParticleCell<Molecule> cell4;
  cell3.addParticle(p7);
  cell3.addParticle(p8);

  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoALoader(cell3, cell3._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoALoader(cell4, cell4._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // This is a basic check for the accumulated values, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.

  // first functors on one cell
  AUTOPAS_OPENMP(parallel) {
    AUTOPAS_OPENMP(sections) {
      AUTOPAS_OPENMP(section)
      functor.SoAFunctorSingle(cell1._particleSoABuffer, newton3);
      AUTOPAS_OPENMP(section)
      functor.SoAFunctorSingle(cell2._particleSoABuffer, newton3);
      AUTOPAS_OPENMP(section)
      functor.SoAFunctorSingle(cell3._particleSoABuffer, newton3);
      AUTOPAS_OPENMP(section)
      functor.SoAFunctorSingle(cell4._particleSoABuffer, newton3);
    }
  }

  // functors on two cells
  AUTOPAS_OPENMP(parallel) {
    AUTOPAS_OPENMP(sections) {
      AUTOPAS_OPENMP(section)
      functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
      AUTOPAS_OPENMP(section)
      functor.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, newton3);
    }
  }
}
