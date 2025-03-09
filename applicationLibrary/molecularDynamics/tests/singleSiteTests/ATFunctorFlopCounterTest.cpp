/**
 * @file ATFunctorFlopCounterTest.cpp
 * @author muehlhaeusser
 * @date 30.07.2024
 */

#include "ATFunctorFlopCounterTest.h"

#include "autopas/AutoPasDecl.h"
#include "autopas/AutoPasImpl.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"
#include "molecularDynamicsLibrary/AxilrodTellerFunctor.h"
#include "testingHelpers/commonTypedefs.h"

// AutoPas is instantiated in AutoPasInstantiations.cpp
// but computeInteractions() versions with countFLOPs == true not,
// so they have to be explicitly instantiated here.
extern template class autopas::AutoPas<Molecule>;
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::AxilrodTellerFunctor<Molecule, false, autopas::FunctorN3Modes::Both, false, /*countFLOPs*/ true> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::AxilrodTellerFunctor<Molecule, false, autopas::FunctorN3Modes::Both, true, /*countFLOPs*/ true> *);

/**
 * Generates a square of four particles, iterates over it with the AxilrodTellerFunctor and checks the values of
 * getNumFLOPs() and getHitRate()
 * @tparam calculateGlobals
 * @param dataLayoutOption
 * @param newton3
 * @param isVerlet
 */
template <bool calculateGlobals>
void ATFunctorFlopCounterTest::testFLOPCounter(autopas::DataLayoutOption dataLayoutOption, bool newton3,
                                               bool isVerlet) {
  const auto isSoA = dataLayoutOption == autopas::DataLayoutOption::soa;

  if (isVerlet or isSoA) {
    autopas::utils::ExceptionHandler::exception(
        "ATFunctorFlopCounterTest::testFLOPCounter Tests are not yet implemented for SoA or VerletSoAs!");
  }

  autopas::AutoPas<Molecule> autoPas;

  autoPas.setBoxMin({0, 0, 0});
  autoPas.setBoxMax({3, 3, 3});
  autoPas.setCutoff(1.1);
  autoPas.setVerletSkinPerTimestep(0.2);
  autoPas.setVerletRebuildFrequency(1);

  // Direct Sum for now to allow newton3 on and off.
  autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
  autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential}, autopas::InteractionTypeOption::triwise);

  if (newton3) {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::enabled}, autopas::InteractionTypeOption::triwise);
  } else {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled}, autopas::InteractionTypeOption::triwise);
  }
  autoPas.setAllowedDataLayouts(std::set<autopas::DataLayoutOption>{dataLayoutOption},
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
    autoPas.addParticle(m, true);
  }

  //  // update container
  //  auto buffer = autoPas.updateContainer();  // buffer is meaningless here

  mdLib::AxilrodTellerFunctor<Molecule, false, autopas::FunctorN3Modes::Both, calculateGlobals, true> atFunctor(
      autoPas.getCutoff());

  autoPas.computeInteractions(&atFunctor);

  // See above for reasoning.
  const auto expectedDistanceCalculations = newton3 ? 4 : 12;

  // See above for reasoning
  const auto expectedNoN3KernelCalls = newton3 ? 0 : 3;
  const auto expectedN3KernelCalls = newton3 ? 1 : 0;

  const auto expectedN3GlobalsCalcs = calculateGlobals ? expectedN3KernelCalls : 0;
  const auto expectedNoN3GlobalsCalcs = calculateGlobals ? expectedNoN3KernelCalls : 0;

  // three distance calculations cost 24 FLOPs, AT kernel calls without Newton3 cost 59 FLOPs, with Newton 3 cost 100
  // FLOPs globals calculations without Newton3 cost 10 FLOPs, with Newton3 cost 24 FLOPs
  constexpr int numFLOPsPerDistanceCalc = 24;
  constexpr int numFLOPsPerNoN3KernelCall = 59;
  constexpr int numFLOPsPerN3KernelCall = 100;
  constexpr int numFLOPsPerNoN3GlobalsCall = 10;
  constexpr int numFLOPsPerN3GlobalsCall = 24;
  const auto expectedFlops =
      expectedDistanceCalculations * numFLOPsPerDistanceCalc + expectedN3KernelCalls * numFLOPsPerN3KernelCall +
      expectedNoN3KernelCalls * numFLOPsPerNoN3KernelCall + expectedN3GlobalsCalcs * numFLOPsPerN3GlobalsCall +
      expectedNoN3GlobalsCalcs * numFLOPsPerNoN3GlobalsCall;
  ASSERT_EQ(expectedFlops, atFunctor.getNumFLOPs());

  const auto expectedHitRate =
      ((double)expectedN3KernelCalls + (double)expectedNoN3KernelCalls) / (double)expectedDistanceCalculations;
  ASSERT_NEAR(expectedHitRate, atFunctor.getHitRate(), 1e-14);
}

template <bool calculateGlobals>
void ATFunctorFlopCounterTest::testFLOPCounterAoSOMP(bool newton3) {
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);
  Molecule p3({0.3, 0.2, 0.1}, {0., 0., 0.}, 2, 0);

  Molecule p4({0., 2., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p5({0.1, 2.2, 0.3}, {0., 0., 0.}, 1, 0);
  Molecule p6({0.3, 2.2, 0.1}, {0., 0., 0.}, 3, 0);

  const double cutoff = 1.;

  mdLib::AxilrodTellerFunctor<Molecule, false, autopas::FunctorN3Modes::Both, calculateGlobals, true> atFunctor(cutoff);

  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    atFunctor.AoSFunctor(p1, p2, p3, newton3);
    AUTOPAS_OPENMP(section)
    atFunctor.AoSFunctor(p4, p5, p6, newton3);
  }
}

/**
 * Tests that the FLOP counts produced are correct by comparing against partially hard-coded values.
 */
TEST_P(ATFunctorFlopCounterTest, testFLOPCountingNoOMP) {
  const auto [dataLayout, newton3, calculateGlobals, isVerlet] = GetParam();
  if (calculateGlobals) {
    ATFunctorFlopCounterTest::testFLOPCounter<true>(dataLayout, newton3, isVerlet);
  } else {
    ATFunctorFlopCounterTest::testFLOPCounter<false>(dataLayout, newton3, isVerlet);
  }
}

/**
 * Tests that FLOP counting has no data races by performing interactions in parallel. With thread sanitizer enabled,
 * this should produce errors. Without thread sanitizer enabled, this test will generally not throw errors.
 */
TEST_P(ATFunctorFlopCounterTest, testFLOPCountingOMP) {
  const auto [dataLayout, newton3, calculateGlobals, isVerlet] = GetParam();
  if (calculateGlobals) {
    if (dataLayout == autopas::DataLayoutOption::aos) {
      ATFunctorFlopCounterTest::testFLOPCounterAoSOMP<true>(newton3);
    } else {
      if (isVerlet) {
        ATFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP<true>(newton3);
      } else {
        ATFunctorFlopCounterTest::testFLOPCounterSoASingleAndPairOMP<true>(newton3);
      }
    }
  } else {
    if (dataLayout == autopas::DataLayoutOption::aos) {
      ATFunctorFlopCounterTest::testFLOPCounterAoSOMP<false>(newton3);
    } else {
      if (isVerlet) {
        ATFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP<false>(newton3);
      } else {
        ATFunctorFlopCounterTest::testFLOPCounterSoASingleAndPairOMP<false>(newton3);
      }
    }
  }
}

/**
 * We test AxilrodTellerFunctor FLOP counting for a combination of data layouts, newton3, calcGlobals, isVerlet.
 *
 * isVerlet is specifically to test the SoA Verlet functor so only relevant with SoA.
 * TODO: Enable SoA and Verlet configurations once they are implemented.
 *
 * @return
 */
INSTANTIATE_TEST_SUITE_P(
    ATFunctorFlopTestSuite, ATFunctorFlopCounterTest,
    /*                               Data Layout                 , n3, calcGlob, isVerlet */
    testing::Values(std::make_tuple(autopas::DataLayoutOption::aos, false, false, false),
                    std::make_tuple(autopas::DataLayoutOption::aos, true, false, false),
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, false, false, false),
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, true, false, false),
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, false, false, true),
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, true, false, true),

                    std::make_tuple(autopas::DataLayoutOption::aos, false, true, false),
                    std::make_tuple(autopas::DataLayoutOption::aos, true, true, false)
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, false, true, false),
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, true, true, false),
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, false, true, true),
                    //                    std::make_tuple(autopas::DataLayoutOption::soa, true, true, true)
                    ));
