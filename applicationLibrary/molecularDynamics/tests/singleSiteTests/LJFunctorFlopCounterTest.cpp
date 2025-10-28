/**
 * @file LJFunctorFlopCounterTest.cpp
 * @author F. Gratl
 * @date 01.06.18
 */

#include "LJFunctorFlopCounterTest.h"

#include "autopas/AutoPasDecl.h"
#include "autopas/AutoPasImpl.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

// AutoPas is instantiated in AutoPasInstantiations.cpp
// but computeInteractions() versions with countFLOPs == true not,
// so they have to be explicitly instantiated here.
extern template class autopas::AutoPas<Molecule>;
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::LJFunctor<Molecule, false, false, autopas::FunctorN3Modes::Both, false, /*countFLOPs*/ true, true> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::LJFunctor<Molecule, false, false, autopas::FunctorN3Modes::Both, true, /*countFLOPs*/ true, true> *);
template bool autopas::AutoPas<Molecule>::computeInteractions(
    mdLib::LJFunctor<Molecule, true, false, autopas::FunctorN3Modes::Both, true, /*countFLOPs*/ true, true> *);

/**
 * Generates a square of four particles, iterates over it with the LJFunctor and checks the values of getNumFLOPs() and
 * getHitRate()
 * @tparam calculateGlobals
 * @tparam applyShift
 * @param dataLayoutOption
 * @param newton3
 * @param isVerlet if true, uses VLC container and traversal. Otherwise used LC.
 */
template <bool calculateGlobals, bool applyShift>
void LJFunctorFlopCounterTest::testFLOPCounter(autopas::DataLayoutOption dataLayoutOption, bool newton3,
                                               bool isVerlet) {
  const auto isSoA = dataLayoutOption == autopas::DataLayoutOption::soa;

  if (isVerlet and not isSoA) {
    autopas::utils::ExceptionHandler::exception(
        "LJFunctorFlopCounterTest::testFLOPCounter using isVerlet without SoA is not allowed!");
  }

  autopas::AutoPas<Molecule> autoPas;

  autoPas.setBoxMin({0, 0, 0});
  autoPas.setBoxMax({3, 3, 3});
  autoPas.setCutoff(1.1);
  autoPas.setVerletSkin(0.2);
  autoPas.setVerletRebuildFrequency(1);
  if (isVerlet) {
    autoPas.setAllowedContainers({autopas::ContainerOption::verletListsCells});
    autoPas.setAllowedTraversals({autopas::TraversalOption::vlc_c18});
  } else {
    autoPas.setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPas.setAllowedTraversals({autopas::TraversalOption::lc_c08});
  }

  if (newton3) {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::enabled});
  } else {
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  }
  autoPas.setAllowedDataLayouts(std::set<autopas::DataLayoutOption>{dataLayoutOption});

  autoPas.init();

  const std::vector<Molecule> molVec{Molecule({0.9, 0.3, 0.9}, {0, 0, 0}, 0), Molecule({0.9, 0.9, 0.9}, {0, 0, 0}, 1),
                                     Molecule({0.9, 1.5, 0.9}, {0, 0, 0}, 2), Molecule({0.1, 2.4, 0.1}, {0, 0, 0}, 3)};

  /**
   * Explanation of molVec choice and resulting numbers of distance calculations and kernel calls.
   *
   * For not isVerlet, we have simple cells: mol 0 & 1 in cell "A", mol 2 & 3 in cell "B"
   *
   * This results in the following interactions:
   * mol 0 <-> 1: SoAFunctorSingle, inside cutoff (dist: 0.6)
   * mol 0 <-> 2: SoAFunctorPair, outside cutoff (dist: 1.2)
   * mol 0 <-> 3: SoAFunctorPair, outside cutoff (dist: ~2.36)
   * mol 1 <-> 2: SoAFunctorPair, inside cutoff (dist: 0.6)
   * mol 1 <-> 3: SoAFunctorPair, outside cutoff (dist: ~1.89)
   * mol 2 <-> 3: SoAFunctorSingle, outside cutoff (dist: ~1.45)
   * This tests all possible combinations of SoAFunctorSingle/Pair with inside/outside cutoff
   *
   * SoAFunctorSingle: 2 distance calls, 1 N3 kernel call (regardless of newton3 enabled/disabled)
   * SoAFunctorPair: N3=enabled : 4 distance calls, 1 N3 kernel call
   *                 N3=disabled: 8 distance calls, 2 NoN3 kernel calls
   *
   * Total: N3=enabled: 6 distance calls, 2 N3 kernel calls
   *        N3=disabled: 10 distance calls, 1 N3 kernel call, 2 NoN3 kernel calls
   *
   * Similarly also tests AoS Functor with inside/outside cutoff:
   * 6 distance calls, 2 kernel calls. Double this with N3=disabled. Kernel calls always match N3=enabled/disabled.
   *
   *
   * For isVerlet-SoA, we have the following neighbor lists (assume N3=disabled):
   * mol 0: 1, 2
   * mol 1: 0, 2
   * mol 2: 0, 1
   * mol 3: Nothing
   * Similar, this tests inside cutoff (0 <-> 1, 1 <-> 2) and outside cutoff (0 <-> 3)
   *
   * We have 3 distance calls, 2 kernel calls. Double this with N3=disabled. Kernel calls always match
   * N3=enabled/disabled.
   */

  for (auto &m : molVec) {
    autoPas.addParticle(m);
  }

  // update container -> build neighbor lists in case of Verlet
  auto buffer = autoPas.updateContainer();  // buffer is meaningless here

  mdLib::LJFunctor<Molecule, applyShift, false, autopas::FunctorN3Modes::Both, calculateGlobals, true, true> ljFunctor(
      autoPas.getCutoff());

  autoPas.computeInteractions(&ljFunctor);

  // See above for reasoning.
  const auto expectedDistanceCalculations =
      isVerlet ? (newton3 ? 3 : 6) : (isSoA ? (newton3 ? 6 : 10) : (newton3 ? 6 : 12));

  // See above for reasoning
  const auto expectedNoN3KernelCalls = newton3 ? 0 : (not isVerlet and isSoA ? 2 : 4);
  const auto expectedN3KernelCalls = newton3 ? 2 : (not isVerlet and isSoA ? 1 : 0);

  const auto expectedN3GlobalsCalcs = calculateGlobals ? expectedN3KernelCalls : 0;
  const auto expectedNoN3GlobalsCalcs = calculateGlobals ? expectedNoN3KernelCalls : 0;

  // distance calculations cost 8 FLOPs, LJ kernel calls without Newton3 cost 15 FLOPs, with Newton 3 cost 18 FLOPs
  // globals calculations without Newton3 cost 8 FLOPs, 9 with shift
  // globals calculations with Newton3 cost 12 FLOPs, 13 with shift (+4 for adding energy and virial to Particle 2)
  constexpr int numFLOPsPerDistanceCalc = 8;
  constexpr int numFLOPsPerNoN3KernelCall = 15;
  constexpr int numFLOPsPerN3KernelCall = 18;
  constexpr int numFLOPsPerNoN3GlobalsCall = applyShift ? 9 : 8;
  constexpr int numFLOPsPerN3GlobalsCall = applyShift ? 13 : 12;
  const auto expectedFlops =
      expectedDistanceCalculations * numFLOPsPerDistanceCalc + expectedN3KernelCalls * numFLOPsPerN3KernelCall +
      expectedNoN3KernelCalls * numFLOPsPerNoN3KernelCall + expectedN3GlobalsCalcs * numFLOPsPerN3GlobalsCall +
      expectedNoN3GlobalsCalcs * numFLOPsPerNoN3GlobalsCall;
  ASSERT_EQ(expectedFlops, ljFunctor.getNumFLOPs());

  const auto expectedHitRate =
      ((double)expectedN3KernelCalls + (double)expectedNoN3KernelCalls) / (double)expectedDistanceCalculations;
  ASSERT_NEAR(expectedHitRate, ljFunctor.getHitRate(), 1e-14);
}

template <bool calculateGlobals, bool applyShift>
void LJFunctorFlopCounterTest::testFLOPCounterAoSOMP(bool newton3) {
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 2., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 2.2, 0.3}, {0., 0., 0.}, 1, 0);

  const double cutoff = 1.;

  mdLib::LJFunctor<Molecule, applyShift, false, autopas::FunctorN3Modes::Both, calculateGlobals, true, true> ljFunctor(
      cutoff);

  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    ljFunctor.AoSFunctor(p1, p2, newton3);
    AUTOPAS_OPENMP(section)
    ljFunctor.AoSFunctor(p3, p4, newton3);
  }
}

template <bool calculateGlobals, bool applyShift>
void LJFunctorFlopCounterTest::testFLOPCounterSoASingleAndPairOMP(bool newton3) {
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p5({1., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p6({1.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p7({1., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p8({1.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  const double cutoff = 1.;

  mdLib::LJFunctor<Molecule, applyShift, false, autopas::FunctorN3Modes::Both, calculateGlobals, true, true> ljFunctor(
      cutoff);

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

  ljFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctor.SoALoader(cell3, cell3._particleSoABuffer, 0, /*skipSoAResize*/ false);
  ljFunctor.SoALoader(cell4, cell4._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // This is a basic check for the accumulated values, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.

  // first functors on one soa
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    ljFunctor.SoAFunctorSingle(cell1._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    ljFunctor.SoAFunctorSingle(cell2._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    ljFunctor.SoAFunctorSingle(cell3._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    ljFunctor.SoAFunctorSingle(cell4._particleSoABuffer, newton3);
  }

  // functors on two soas
  AUTOPAS_OPENMP(parallel sections) {
    AUTOPAS_OPENMP(section)
    ljFunctor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    AUTOPAS_OPENMP(section)
    ljFunctor.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, newton3);
  }
}

template <bool calculateGlobals, bool applyShift>
void LJFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP(bool newton3) {
  Molecule p0({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p1({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p2({0., 2., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p3({0.1, 2.2, 0.3}, {0., 0., 0.}, 1, 0);

  // generate neighbor lists
  std::array<std::vector<size_t, autopas::AlignedAllocator<size_t>>, 4> neighborLists;
  neighborLists[0].push_back(1);  // p0 has neighbor p1
  neighborLists[1].push_back(0);  // p1 has neighbor p0
  neighborLists[2].push_back(3);  // p2 has neighbor p3
  neighborLists[3].push_back(2);  // p3 has neighbor p2

  const double cutoff = 1.;

  mdLib::LJFunctor<Molecule, applyShift, false, autopas::FunctorN3Modes::Both, calculateGlobals, true, true> ljFunctor(
      cutoff);

  autopas::FullParticleCell<Molecule> cell;
  cell.addParticle(p0);
  cell.addParticle(p1);
  cell.addParticle(p2);
  cell.addParticle(p3);

  ljFunctor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);

  // This is a basic check for the accumulated values, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.

  AUTOPAS_OPENMP(parallel) {
    AUTOPAS_OPENMP(sections) {
      AUTOPAS_OPENMP(section)
      ljFunctor.SoAFunctorVerlet(cell._particleSoABuffer, 0, neighborLists[0], newton3);
      AUTOPAS_OPENMP(section)
      ljFunctor.SoAFunctorVerlet(cell._particleSoABuffer, 2, neighborLists[2], newton3);
    }
  }
}

/**
 * Tests that the FLOP counts produced are correct by comparing against partially hard-coded values.
 */
TEST_P(LJFunctorFlopCounterTest, testFLOPCountingNoOMP) {
  const auto [dataLayout, newton3, calculateGlobals, applyShift, isVerlet] = GetParam();
  if (calculateGlobals and applyShift) {
    LJFunctorFlopCounterTest::testFLOPCounter<true, true>(dataLayout, newton3, isVerlet);
  } else if (calculateGlobals and not applyShift) {
    LJFunctorFlopCounterTest::testFLOPCounter<true, false>(dataLayout, newton3, isVerlet);
  } else {
    LJFunctorFlopCounterTest::testFLOPCounter<false, false>(dataLayout, newton3, isVerlet);
  }
}

/**
 * Tests that FLOP counting has no data races by performing interactions in parallel. With thread sanitizer enabled,
 * this should produce errors. Without thread sanitizer enabled, this test will generally not throw errors.
 */
TEST_P(LJFunctorFlopCounterTest, testFLOPCountingOMP) {
  const auto [dataLayout, newton3, calculateGlobals, applyShift, isVerlet] = GetParam();
  if (calculateGlobals and applyShift) {
    if (dataLayout == autopas::DataLayoutOption::aos) {
      LJFunctorFlopCounterTest::testFLOPCounterAoSOMP<true, true>(newton3);
    } else {
      if (isVerlet) {
        LJFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP<true, true>(newton3);
      } else {
        LJFunctorFlopCounterTest::testFLOPCounterSoASingleAndPairOMP<true, true>(newton3);
      }
    }
  } else if (calculateGlobals and not applyShift) {
    if (dataLayout == autopas::DataLayoutOption::aos) {
      LJFunctorFlopCounterTest::testFLOPCounterAoSOMP<true, false>(newton3);
    } else {
      if (isVerlet) {
        LJFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP<true, false>(newton3);
      } else {
        LJFunctorFlopCounterTest::testFLOPCounterSoASingleAndPairOMP<true, false>(newton3);
      }
    }
  } else {
    if (dataLayout == autopas::DataLayoutOption::aos) {
      LJFunctorFlopCounterTest::testFLOPCounterAoSOMP<false, false>(newton3);
    } else {
      if (isVerlet) {
        LJFunctorFlopCounterTest::testFLOPCounterSoAVerletOMP<false, false>(newton3);
      } else {
        LJFunctorFlopCounterTest::testFLOPCounterSoASingleAndPairOMP<false, false>(newton3);
      }
    }
  }
}

/**
 * We test LJFunctor FLOP counting for a combination of data layouts, newton3, calcGlobals, applyShift, isVerlet.
 *
 * applyShift is only relevant when calculateGlobals=true.
 * isVerlet is specifically to test the SoA Verlet functor so only relevant with SoA.
 *
 * @return
 */
INSTANTIATE_TEST_SUITE_P(
    LJFunctorFlopTestSuite, LJFunctorFlopCounterTest,
    /*                               Data Layout                 , n3, calcGlob, appShift, isVerlet */
    testing::Values(std::make_tuple(autopas::DataLayoutOption::aos, false, false, false, false),
                    std::make_tuple(autopas::DataLayoutOption::aos, true, false, false, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, false, false, false, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, true, false, false, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, false, false, false, true),
                    std::make_tuple(autopas::DataLayoutOption::soa, true, false, false, true),

                    std::make_tuple(autopas::DataLayoutOption::aos, false, true, false, false),
                    std::make_tuple(autopas::DataLayoutOption::aos, true, true, false, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, false, true, false, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, true, true, false, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, false, true, false, true),
                    std::make_tuple(autopas::DataLayoutOption::soa, true, true, false, true),

                    std::make_tuple(autopas::DataLayoutOption::aos, false, true, true, false),
                    std::make_tuple(autopas::DataLayoutOption::aos, true, true, true, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, false, true, true, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, true, true, true, false),
                    std::make_tuple(autopas::DataLayoutOption::soa, false, true, true, true),
                    std::make_tuple(autopas::DataLayoutOption::soa, true, true, true, true)));
