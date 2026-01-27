/**
 * @file ATMFunctorFlopCounterTest.cpp
 * @author muehlhaeusser
 * @date 30.07.2024
 */

#include "ATMFunctorFlopCounterTest.h"

#include "autopas/AutoPasDecl.h"
#include "autopas/AutoPasImpl.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"
#include "molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h"
#include "testingHelpers/commonTypedefs.h"

// AutoPas is instantiated in AutoPasInstantiations.cpp
// but computeInteractions() versions with countFLOPs == true not,
// so they have to be explicitly instantiated here.
extern template class autopas::AutoPas<Molecule>;
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, false, /*countFLOPs*/ true> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true, /*countFLOPs*/ true> *);

/**
 * Generates a square of four particles, iterates over it with the AxilrodTellerMutoFunctor and checks the values of
 * getNumFLOPs() and getHitRate()
 * @tparam calculateGlobals
 * @param newton3
 */
template <bool calculateGlobals>
void ATMFunctorFlopCounterTest::testFLOPCounterAoS(bool newton3) {
  autopas::AutoPas<Molecule> autoPas;

  autoPas.setBoxMin({0, 0, 0});
  autoPas.setBoxMax({3, 3, 3});
  autoPas.setCutoff(1.1);
  autoPas.setVerletSkin(0.2);
  autoPas.setVerletRebuildFrequency(1);

  // Direct Sum for now to allow newton3 on and off.
  autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
  autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential}, autopas::InteractionTypeOption::triwise);

  if (newton3) {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::enabled}, autopas::InteractionTypeOption::triwise);
  } else {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled}, autopas::InteractionTypeOption::triwise);
  }
  autoPas.setAllowedDataLayouts(std::set<autopas::DataLayoutOption>{autopas::DataLayoutOption::aos},
                                autopas::InteractionTypeOption::triwise);

  autoPas.setAllowedInteractionTypeOptions({autopas::InteractionTypeOption::triwise});
  autoPas.init();

  const std::vector<Molecule> molVec{Molecule({0.2, 0.2, 0.2}, {0, 0, 0}, 0), Molecule({1.0, 0.2, 0.2}, {0, 0, 0}, 1),
                                     Molecule({1.0, 0.8, 0.2}, {0, 0, 0}, 2), Molecule({0.2, 0.2, 2.5}, {0, 0, 0}, 3)};

  /**
   * Explanation of molVec choice and resulting numbers of distance calculations and kernel calls.
   *
   * Currently, this test is only meant to test flop counting for the AoSFunctor
   *
   * Interactions within the cutoff: {0, 1, 2}
   * Interactions outside the cutoff: {0, 1, 3}, {0, 2, 3}, {1, 2, 3}
   *
   * Distance Calls: 4 with N3, 12 without N3
   * Kernel Calls: 1 with N3, 3 without N3
   */

  for (auto &m : molVec) {
    autoPas.addParticle(m);
  }

  //  update container
  //  auto buffer = autoPas.updateContainer();  // buffer is meaningless here

  mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, calculateGlobals, true> atmFunctor(
      autoPas.getCutoff());

  autoPas.computeInteractions(&atmFunctor);

  // See above for reasoning.
  const auto expectedDistanceCalculations = newton3 ? 4 : 12;

  // See above for reasoning
  const auto expectedNoN3KernelCalls = newton3 ? 0 : 3;
  const auto expectedN3KernelCalls = newton3 ? 1 : 0;

  const auto expectedN3GlobalsCalcs = calculateGlobals ? expectedN3KernelCalls : 0;
  const auto expectedNoN3GlobalsCalcs = calculateGlobals ? expectedNoN3KernelCalls : 0;

  // three distance calculations cost 24 FLOPs, AT kernel calls without Newton3 cost 59 FLOPs, with Newton 3 cost 100
  // FLOPs globals calculations without Newton3 cost 13 FLOPs, with Newton3 cost 33 FLOPs
  constexpr int numFLOPsPerDistanceCalc = 24;
  constexpr int numFLOPsPerNoN3KernelCall = 59;
  constexpr int numFLOPsPerN3KernelCall = 100;
  constexpr int numFLOPsPerNoN3GlobalsCall = 13;
  constexpr int numFLOPsPerN3GlobalsCall = 33;
  const auto expectedFlops =
      expectedDistanceCalculations * numFLOPsPerDistanceCalc + expectedN3KernelCalls * numFLOPsPerN3KernelCall +
      expectedNoN3KernelCalls * numFLOPsPerNoN3KernelCall + expectedN3GlobalsCalcs * numFLOPsPerN3GlobalsCall +
      expectedNoN3GlobalsCalcs * numFLOPsPerNoN3GlobalsCall;
  ASSERT_EQ(expectedFlops, atmFunctor.getNumFLOPs());

  const auto expectedHitRate =
      ((double)expectedN3KernelCalls + (double)expectedNoN3KernelCalls) / (double)expectedDistanceCalculations;
  ASSERT_NEAR(expectedHitRate, atmFunctor.getHitRate(), 1e-14);
}

template <bool calculateGlobals>
void ATMFunctorFlopCounterTest::testFLOPCounterAoSOMP(bool newton3) {
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);
  Molecule p3({0.3, 0.2, 0.1}, {0., 0., 0.}, 2, 0);

  Molecule p4({0., 2., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p5({0.1, 2.2, 0.3}, {0., 0., 0.}, 1, 0);
  Molecule p6({0.3, 2.2, 0.1}, {0., 0., 0.}, 3, 0);

  const double cutoff = 1.;

  mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, calculateGlobals, true> atmFunctor(
      cutoff);

  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    atmFunctor.AoSFunctor(p1, p2, p3, newton3);
    AUTOPAS_OPENMP(section)
    atmFunctor.AoSFunctor(p4, p5, p6, newton3);
  }
}

/**
 * Generates four particles and places them into either one, two or three cells.
 * Then iterates over it with the AxilrodTellerMutoFunctor and checks the values of getNumFLOPs() and getHitRate()
 * @tparam calculateGlobals
 * @param newton3
 * @param functorFunctionType
 */
template <bool calculateGlobals>
void ATMFunctorFlopCounterTest::testFLOPCounterSoA(bool newton3, FunctorFunction functorFunctionType) {
  if (functorFunctionType == FunctorFunction::verlet) {
    autopas::utils::ExceptionHandler::exception(
        "ATMFunctorFlopCounterTest::testFLOPCounterAoS Tests are not yet implemented for VerletSoAs!");
  }

  constexpr double cutoff = 1.1;

  // const auto particlePositions = getParticlePositions(functorFunctionType);
  constexpr std::array p0 = {0.2, 0.2, 0.2};
  constexpr std::array p1 = {1.0, 0.2, 0.2};
  constexpr std::array p2 = {1.0, 0.8, 0.2};
  constexpr std::array p3 = {0.2, 0.2, 2.5};
  const std::vector molVec{Molecule(p0, {0, 0, 0}, 0), Molecule(p1, {0, 0, 0}, 1), Molecule(p2, {0, 0, 0}, 2),
                           Molecule(p3, {0, 0, 0}, 3)};

  size_t expectedNoN3KernelCalls{}, expectedN3KernelCalls{}, expectedN3GlobalsCalcs{}, expectedNoN3GlobalsCalcs{};
  size_t maxNumDistanceCalculations = 12;
  size_t expectedNumTriplets = 4;

  mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, calculateGlobals, true> atmFunctor(
      cutoff);
  atmFunctor.initTraversal();
  switch (functorFunctionType) {
    case soaSingle: {
      FMCell cell0;
      for (auto &m : molVec) {
        cell0.addParticle(m);
      }
      atmFunctor.SoALoader(cell0, cell0._particleSoABuffer, 0, false);
      atmFunctor.SoAFunctorSingle(cell0._particleSoABuffer, newton3);
      expectedN3KernelCalls = 1;
      expectedN3GlobalsCalcs = calculateGlobals ? 1 : 0;
      break;
    }
    case soaPair: {
      FMCell cell0, cell1;
      cell0.addParticle(molVec[0]);
      cell0.addParticle(molVec[1]);
      cell1.addParticle(molVec[2]);
      cell1.addParticle(molVec[3]);

      atmFunctor.SoALoader(cell0, cell0._particleSoABuffer, 0, false);
      atmFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
      atmFunctor.SoAFunctorPair(cell0._particleSoABuffer, cell1._particleSoABuffer, newton3);

      if (not newton3) {
        atmFunctor.SoAFunctorPair(cell1._particleSoABuffer, cell0._particleSoABuffer, newton3);
        expectedNumTriplets = 8;
        maxNumDistanceCalculations = 24;
        expectedNoN3KernelCalls = 1;
        expectedNoN3GlobalsCalcs = calculateGlobals ? 1 : 0;
        expectedN3KernelCalls = 1;
        expectedN3GlobalsCalcs = calculateGlobals ? 1 : 0;
      } else {
        maxNumDistanceCalculations = 12;
        expectedN3KernelCalls = 1;
        expectedN3GlobalsCalcs = calculateGlobals ? 1 : 0;
      }
      break;
    }
    case soaTriple: {
      FMCell cell0, cell1, cell2;
      // Distribute the valid triplet over 3 cells {0, 1, 2}
      cell0.addParticle(molVec[0]);
      cell1.addParticle(molVec[1]);
      cell2.addParticle(molVec[2]);
      cell0.addParticle(molVec[3]);
      atmFunctor.SoALoader(cell0, cell0._particleSoABuffer, 0, false);
      atmFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
      atmFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
      atmFunctor.SoAFunctorTriple(cell0._particleSoABuffer, cell1._particleSoABuffer, cell2._particleSoABuffer,
                                  newton3);

      if (not newton3) {
        atmFunctor.SoAFunctorTriple(cell1._particleSoABuffer, cell0._particleSoABuffer, cell2._particleSoABuffer,
                                    newton3);
        atmFunctor.SoAFunctorTriple(cell2._particleSoABuffer, cell0._particleSoABuffer, cell1._particleSoABuffer,
                                    newton3);
        expectedNumTriplets = 6;
        maxNumDistanceCalculations = 36;
        expectedNoN3KernelCalls = 3;
        expectedNoN3GlobalsCalcs = calculateGlobals ? 3 : 0;
      } else {
        expectedNumTriplets = 2;
        expectedN3KernelCalls = 1;
        expectedN3GlobalsCalcs = calculateGlobals ? 1 : 0;
      }
      break;
    }
    default:
      break;
  }

  atmFunctor.endTraversal(newton3);

  const auto countedFlops = atmFunctor.getNumFLOPs();
  const auto countedHitRate = atmFunctor.getHitRate();

  // One distance calculations costs 8 FLOPs, AT kernel call without Newton3 costs 59 FLOPs, with Newton 3 costs 100
  // FLOPs globals calculations without Newton3 cost 13 FLOPs, with Newton3 cost 33 FLOPs
  constexpr int numFLOPsPerDistanceCalc = 8;
  constexpr int numFLOPsPerNoN3KernelCall = 59;
  constexpr int numFLOPsPerN3KernelCall = 100;
  constexpr int numFLOPsPerNoN3GlobalsCall = 13;
  constexpr int numFLOPsPerN3GlobalsCall = 33;
  // Maximum number of FLOPs computed with the maximum number of distance calculations
  const auto expectedFlopsMax =
      maxNumDistanceCalculations * numFLOPsPerDistanceCalc + expectedN3KernelCalls * numFLOPsPerN3KernelCall +
      expectedNoN3KernelCalls * numFLOPsPerNoN3KernelCall + expectedN3GlobalsCalcs * numFLOPsPerN3GlobalsCall +
      expectedNoN3GlobalsCalcs * numFLOPsPerNoN3GlobalsCall;
  // Minimum number of FLOPs assumes just one distance calculation per triplet
  const auto expectedFlopsMin =
      expectedNumTriplets * numFLOPsPerDistanceCalc + expectedN3KernelCalls * numFLOPsPerN3KernelCall +
      expectedNoN3KernelCalls * numFLOPsPerNoN3KernelCall + expectedN3GlobalsCalcs * numFLOPsPerN3GlobalsCall +
      expectedNoN3GlobalsCalcs * numFLOPsPerNoN3GlobalsCall;

  // We cannot predict the exact number of FLOPs as distance checks can be 1, 2 or 3 per triplet
  ASSERT_LE(countedFlops, expectedFlopsMax) << "Tested for " << to_string(functorFunctionType);
  ASSERT_GE(countedFlops, expectedFlopsMin) << "Tested for " << to_string(functorFunctionType);

  const auto expectedHitRate =
      (static_cast<double>(expectedN3KernelCalls) + static_cast<double>(expectedNoN3KernelCalls)) /
      static_cast<double>(expectedNumTriplets);
  ASSERT_NEAR(countedHitRate, expectedHitRate, 1e-14) << "Tested for " << to_string(functorFunctionType);
}

template <bool calculateGlobals>
void ATMFunctorFlopCounterTest::testFLOPCounterSoAOMP(bool newton3) {
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p5({1., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p6({1.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p7({1., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p8({1.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p9({1., 0., 1.0}, {0., 0., 0.}, 0, 0);
  Molecule p10({1.1, 0.2, 1.3}, {0., 0., 0.}, 1, 0);

  Molecule p11({1., 1., 1.}, {0., 0., 0.}, 0, 0);
  Molecule p12({1.1, 1.2, 1.3}, {0., 0., 0.}, 1, 0);

  const double cutoff = 1.;

  mdLib::AxilrodTellerMutoFunctor<Molecule, false, autopas::FunctorN3Modes::Both, calculateGlobals, true> atmFunctor(
      cutoff);

  autopas::FullParticleCell<Molecule> cell1;
  cell1.addParticle(p1);
  cell1.addParticle(p2);
  cell1.addParticle(p3);

  autopas::FullParticleCell<Molecule> cell2;
  cell2.addParticle(p4);
  cell2.addParticle(p5);

  autopas::FullParticleCell<Molecule> cell3;
  cell3.addParticle(p6);

  autopas::FullParticleCell<Molecule> cell4;
  cell3.addParticle(p7);
  cell3.addParticle(p8);
  cell3.addParticle(p9);

  autopas::FullParticleCell<Molecule> cell5;
  cell4.addParticle(p10);
  cell4.addParticle(p11);

  autopas::FullParticleCell<Molecule> cell6;
  cell5.addParticle(p12);

  atmFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctor.SoALoader(cell3, cell3._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctor.SoALoader(cell4, cell4._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctor.SoALoader(cell5, cell5._particleSoABuffer, 0, /*skipSoAResize*/ false);
  atmFunctor.SoALoader(cell6, cell6._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // This is a basic check for the accumulated values, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.

  // SoA functor for a single cell (only cells with 3+ particles relevant)
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    atmFunctor.SoAFunctorSingle(cell1._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    atmFunctor.SoAFunctorSingle(cell4._particleSoABuffer, newton3);
  }

  // SoA functors on two cells
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    atmFunctor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    atmFunctor.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    atmFunctor.SoAFunctorPair(cell5._particleSoABuffer, cell6._particleSoABuffer, newton3);
  }

  // SoA functors on three cells
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    atmFunctor.SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    atmFunctor.SoAFunctorTriple(cell4._particleSoABuffer, cell5._particleSoABuffer, cell6._particleSoABuffer, newton3);
  }
}

/**
 * Tests that the FLOP counts produced are correct by comparing against partially hard-coded values.
 */
TEST_P(ATMFunctorFlopCounterTest, testFLOPCountingAoSNoOMP) {
  const auto [newton3, calculateGlobals] = GetParam();
  if (calculateGlobals) {
    ATMFunctorFlopCounterTest::testFLOPCounterAoS<true>(newton3);
  } else {
    ATMFunctorFlopCounterTest::testFLOPCounterAoS<false>(newton3);
  }
}

/**
 * Tests that the FLOP counts produced are correct by comparing against partially hard-coded values.
 */
TEST_P(ATMFunctorFlopCounterTest, testFLOPCountingSoANoOMP) {
  const auto [newton3, calculateGlobals] = GetParam();
  auto functorFunctionOptions = {ATMFunctorFlopCounterTest::FunctorFunction::soaSingle,
                                 ATMFunctorFlopCounterTest::FunctorFunction::soaPair,
                                 ATMFunctorFlopCounterTest::FunctorFunction::soaTriple};

  for (auto functorOption : functorFunctionOptions) {
    if (calculateGlobals) {
      ATMFunctorFlopCounterTest::testFLOPCounterSoA<true>(newton3, functorOption);
    } else {
      ATMFunctorFlopCounterTest::testFLOPCounterSoA<false>(newton3, functorOption);
    }
  }
}

/**
 * Tests that FLOP counting has no data races by performing interactions in parallel. With thread sanitizer enabled,
 * this should produce errors. Without thread sanitizer enabled, this test will generally not throw errors.
 */
TEST_P(ATMFunctorFlopCounterTest, testFLOPCountingAoSOMP) {
  const auto [newton3, calculateGlobals] = GetParam();
  if (calculateGlobals) {
    ATMFunctorFlopCounterTest::testFLOPCounterAoSOMP<true>(newton3);
  } else {
    ATMFunctorFlopCounterTest::testFLOPCounterAoSOMP<false>(newton3);
  }
}

/**
 * Tests that FLOP counting has no data races by performing interactions in parallel. With thread sanitizer enabled,
 * this should produce errors. Without thread sanitizer enabled, this test will generally not throw errors.
 */
TEST_P(ATMFunctorFlopCounterTest, testFLOPCountingSoAOMP) {
  const auto [newton3, calculateGlobals] = GetParam();
  auto functorFunctionOptions = {ATMFunctorFlopCounterTest::FunctorFunction::soaSingle,
                                 ATMFunctorFlopCounterTest::FunctorFunction::soaPair,
                                 ATMFunctorFlopCounterTest::FunctorFunction::soaTriple};

  for (auto functorOption : functorFunctionOptions) {
    if (calculateGlobals) {
      if (functorOption == verlet) {
        ATMFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP<true>(newton3);
      } else {
        ATMFunctorFlopCounterTest::testFLOPCounterSoAOMP<true>(newton3);
      }
    } else {
      if (functorOption == verlet) {
        ATMFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP<false>(newton3);
      } else {
        ATMFunctorFlopCounterTest::testFLOPCounterSoAOMP<false>(newton3);
      }
    }
  }
}

/**
 * We test AxilrodTellerMutoFunctor FLOP counting for a combination newton3 and calcGlobals.
 *
 * TODO: Enable Verlet configurations once they are implemented.
 */
INSTANTIATE_TEST_SUITE_P(ATMFunctorFlopTestSuite, ATMFunctorFlopCounterTest,
                         /*  n3, calcGlob */
                         testing::Values(std::make_tuple(false, false), std::make_tuple(true, false),
                                         std::make_tuple(false, true), std::make_tuple(true, true)));
