/**
 * @file CellFunctorTest3B.cpp
 * @author muehlhaeusser
 * @date 29.08.23
 */

#include "CellFunctorTest3B.h"

#include <sstream>

#include "molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h"

// Type aliases via inheritance for more readable test names (using declarations do not work for this)
struct CellFunctor_AoS_NoN3_NoBi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, false> {
  CellFunctor_AoS_NoN3_NoBi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::aos, false) {}
};
struct CellFunctor_AoS_NoN3_Bi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, true> {
  CellFunctor_AoS_NoN3_Bi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::aos, false) {}
};
struct CellFunctor_AoS_N3_NoBi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, false> {
  CellFunctor_AoS_N3_NoBi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::aos, true) {}
};
struct CellFunctor_AoS_N3_Bi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, true> {
  CellFunctor_AoS_N3_Bi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::aos, true) {}
};

struct CellFunctor_SoA_NoN3_NoBi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, false> {
  CellFunctor_SoA_NoN3_NoBi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::soa, false) {}
};
struct CellFunctor_SoA_NoN3_Bi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, true> {
  CellFunctor_SoA_NoN3_Bi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::soa, false) {}
};
struct CellFunctor_SoA_N3_NoBi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, false> {
  CellFunctor_SoA_N3_NoBi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::soa, true) {}
};
struct CellFunctor_SoA_N3_Bi
    : public autopas::internal::CellFunctor3B<autopas::FullParticleCell<Molecule>,
                                              mdLib::AxilrodTellerMutoFunctor<Molecule>, true> {
  CellFunctor_SoA_N3_Bi(mdLib::AxilrodTellerMutoFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor3B(f, sortingCutoff, autopas::DataLayoutOption::soa, true) {}
};

/**
 * All relevant CellFunctor configurations which are used to instantiate the typed test cases
 * testOwnedAndHaloCellInteractionTriple, testOwnedAndHaloCellInteractionPair and testOwnedAndHaloCellInteractionSingle
 *
 */
// clang-format off
using CellFTestingTypes = ::testing::Types<CellFunctor_AoS_NoN3_NoBi,
                                           CellFunctor_AoS_NoN3_Bi,
                                           CellFunctor_AoS_N3_NoBi,
                                           CellFunctor_AoS_N3_Bi>;
                                           //CellFunctor_SoA_NoN3_NoBi,
                                           //CellFunctor_SoA_NoN3_Bi,
                                           //CellFunctor_SoA_N3_NoBi,
                                           //CellFunctor_SoA_N3_Bi>;
// clang-format on

/**
 * Helper for testOwnedAndHaloCellInteractionSingle that creates particles and
 * cells with given OwnershipStates, executes CellFunctor::processCellPair or CellFunctor::processCell and returns the
 * calculated forces on particles.
 *
 * @tparam T The type of CellFunctor3B
 * @param cellFunctor
 * @param ownershipParticle1 OwnershipState for particle 1
 * @param ownershipParticle2 OwnershipState for particle 2
 * @param ownershipParticle3 OwnershipState for particle 3
 * @param ownershipCell OwnershipState of the cell
 * @param dataLayout AoS or SoA
 * @param ATMFunctor The functor to evaluate forces
 * @return The force magnitudes exerted on the 3 particles <Force[p1], Force[p2], Force[p3]>
 */
template <typename T>
std::tuple<double, double, double> doSingleCellInteraction(
    T &cellFunctor, const autopas::OwnershipState ownershipParticle1, const autopas::OwnershipState ownershipParticle2,
    const autopas::OwnershipState ownershipParticle3, const autopas::OwnershipState ownershipCell,
    const autopas::DataLayoutOption dataLayout, mdLib::AxilrodTellerMutoFunctor<Molecule> &ATMFunctor) {
  using autopas::utils::ArrayMath::L2Norm;
  //   create three particles
  Molecule p1({0.8, 0.5, 0.5}, {0., 0., 0.}, 0);
  p1.setOwnershipState(ownershipParticle1);
  Molecule p2({1.2, 0.5, 0.5}, {0., 0., 0.}, 1);
  p2.setOwnershipState(ownershipParticle2);
  Molecule p3({1.0, 1.0, 0.5}, {0., 0., 0.}, 2);
  p3.setOwnershipState(ownershipParticle3);

  // create cell
  autopas::FullParticleCell<Molecule> cell1({1., 1., 1.});

  cell1.setPossibleParticleOwnerships(ownershipCell);
  cell1.addParticle(p1);
  cell1.addParticle(p2);
  cell1.addParticle(p3);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ATMFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }

  cellFunctor.processCell(cell1);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ATMFunctor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
  }

  return std::tuple<double, double, double>{L2Norm(cell1[0].getF()), L2Norm(cell1[1].getF()), L2Norm(cell1[2].getF())};
}

/**
 * Helper for testOwnedAndHaloCellInteractionPair that creates three particles and two cells with given OwnershipStates,
 * with one cell containing 1 particle and the other one 2 particles. Then executes CellFunctor3B::processCellPair and
 * returns the calculated forces on particles.
 *
 * @tparam T The type of CellFunctor3B
 * @param cellFunctor
 * @param ownershipParticle1 OwnershipState for particle 1
 * @param ownershipParticle2 OwnershipState for particle 2
 * @param ownershipParticle3 OwnershipState for particle 3
 * @param ownershipCell1 OwnershipState of the cell 1
 * @param ownershipCell2 OwnershipState of the cell 2
 * @param dataLayout AoS or SoA
 * @param ATMFunctor The functor to evaluate forces
 * @param particle2InCell1 Boolean if the second particle is in cell 1
 * @return The force magnitudes exerted on the 3 particles <Force[p1], Force[p2], Force[p3]>
 */
template <typename T>
std::tuple<double, double, double> doPairCellInteraction(
    T &cellFunctor, const autopas::OwnershipState ownershipParticle1, const autopas::OwnershipState ownershipParticle2,
    const autopas::OwnershipState ownershipParticle3, const autopas::OwnershipState ownershipCell1,
    const autopas::OwnershipState ownershipCell2, const autopas::DataLayoutOption dataLayout,
    mdLib::AxilrodTellerMutoFunctor<Molecule> &ATMFunctor, const bool particle2InCell1) {
  using autopas::utils::ArrayMath::L2Norm;
  //   create three particles
  Molecule p1({0.8, 0.5, 0.5}, {0., 0., 0.}, 0);
  p1.setOwnershipState(ownershipParticle1);
  Molecule p2({1.2, 0.5, 0.5}, {0., 0., 0.}, 1);
  p2.setOwnershipState(ownershipParticle2);
  Molecule p3({1.0, 1.0, 0.5}, {0., 0., 0.}, 2);
  p3.setOwnershipState(ownershipParticle3);

  // create three cells
  autopas::FullParticleCell<Molecule> cell1({1., 1., 1.});
  autopas::FullParticleCell<Molecule> cell2({1., 1., 1.});

  cell1.setPossibleParticleOwnerships(ownershipCell1);
  cell1.addParticle(p1);
  cell2.setPossibleParticleOwnerships(ownershipCell2);
  if (particle2InCell1) {
    cell1.addParticle(p2);
  } else {
    cell2.addParticle(p2);
  }
  cell2.addParticle(p3);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ATMFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ATMFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }

  cellFunctor.processCellPair(cell1, cell2);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ATMFunctor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    ATMFunctor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
  }

  if (particle2InCell1) {
    return std::tuple<double, double, double>{L2Norm(cell1[0].getF()), L2Norm(cell1[1].getF()),
                                              L2Norm(cell2[0].getF())};
  } else {
    return std::tuple<double, double, double>{L2Norm(cell1[0].getF()), L2Norm(cell2[0].getF()),
                                              L2Norm(cell2[1].getF())};
  }
}

/**
 * Helper for testOwnedAndHaloCellInteractionTriple that creates three particles and three cells with given
 * OwnershipStates, with one particle per cell. Then executes CellFunctor3B::processCellTriple and returns the
 * calculated forces on particles.
 *
 * @tparam T The type of CellFunctor3B
 * @param cellFunctor
 * @param ownershipParticle1 OwnershipState for particle 1
 * @param ownershipParticle2 OwnershipState for particle 2
 * @param ownershipParticle3 OwnershipState for particle 3
 * @param ownershipCell1 OwnershipState of the cell 1
 * @param ownershipCell2 OwnershipState of the cell 2
 * @param ownershipCell3 OwnershipState of the cell 3
 * @param dataLayout AoS or SoA
 * @param ATMFunctor The functor to evaluate forces
 * @return The force magnitudes exerted on the 3 particles <Force[p1], Force[p2], Force[p3]>
 */
template <typename T>
std::tuple<double, double, double> doTripleCellInteraction(
    T &cellFunctor, const autopas::OwnershipState ownershipParticle1, const autopas::OwnershipState ownershipParticle2,
    const autopas::OwnershipState ownershipParticle3, const autopas::OwnershipState ownershipCell1,
    const autopas::OwnershipState ownershipCell2, const autopas::OwnershipState ownershipCell3,
    const autopas::DataLayoutOption dataLayout, mdLib::AxilrodTellerMutoFunctor<Molecule> &ATMFunctor) {
  using autopas::utils::ArrayMath::L2Norm;
  //   create three particles
  Molecule p1({0.8, 0.5, 0.5}, {0., 0., 0.}, 0);
  p1.setOwnershipState(ownershipParticle1);
  Molecule p2({1.2, 0.5, 0.5}, {0., 0., 0.}, 1);
  p2.setOwnershipState(ownershipParticle2);
  Molecule p3({1.0, 1.0, 0.5}, {0., 0., 0.}, 2);
  p3.setOwnershipState(ownershipParticle3);

  // create three cells
  autopas::FullParticleCell<Molecule> cell1({1., 1., 1.});
  autopas::FullParticleCell<Molecule> cell2({1., 1., 1.});
  autopas::FullParticleCell<Molecule> cell3({1., 1., 1.});

  cell1.setPossibleParticleOwnerships(ownershipCell1);
  cell1.addParticle(p1);
  cell2.setPossibleParticleOwnerships(ownershipCell2);
  cell2.addParticle(p2);
  cell3.setPossibleParticleOwnerships(ownershipCell3);
  cell3.addParticle(p3);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ATMFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ATMFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
    ATMFunctor.SoALoader(cell3, cell3._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }

  cellFunctor.processCellTriple(cell1, cell2, cell3);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ATMFunctor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    ATMFunctor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    ATMFunctor.SoAExtractor(cell3, cell3._particleSoABuffer, 0);
  }

  return std::tuple<double, double, double>{L2Norm(cell1[0].getF()), L2Norm(cell2[0].getF()), L2Norm(cell3[0].getF())};
}

TYPED_TEST_SUITE_P(CellFunctorTest3B);

/**
 * Tests if force calculation is skipped in the CellFunctor3B with a pure halo cell. Checks if force calculations are
 * done for a cell that can contain owned particles.
 */
TYPED_TEST_P(CellFunctorTest3B, testOwnedAndHaloCellInteractionSingle) {
  using CellFunctorType = TypeParam;

  // shorthands for readability
  constexpr autopas::OwnershipState owned = autopas::OwnershipState::owned;
  constexpr autopas::OwnershipState halo = autopas::OwnershipState::halo;
  constexpr autopas::OwnershipState ownedOrHalo = autopas::OwnershipState::owned | autopas::OwnershipState::halo;

  // helper function to create a summarizing error string
  const auto createErrorString = [&](const auto ownershipCell, const auto ownershipP1, const auto ownershipP2,
                                     const auto ownershipP3, const int numParticle, const bool actualForceZero) {
    std::stringstream errorStr;
    errorStr << "Particle " << numParticle << (actualForceZero ? " does not experience" : " experiences") << " force."
             << "\nCell1: " << ownershipCell << "\nParticle 1 in Cell1: " << ownershipP1
             << "\nParticle 2 in Cell1: " << ownershipP2 << "\nParticle 3 in Cell1: " << ownershipP3;
    return errorStr.str();
  };

  // Test with and without sorting
  for (const auto sortingThreshold : {0, 100}) {
    // Test all reasonable combinations of owned / halo particles and cells
    for (const auto ownershipParticle1 : {owned, halo}) {
      for (const auto ownershipParticle2 : {owned, halo}) {
        for (const auto ownershipParticle3 : {owned, halo}) {
          for (const auto ownershipCell1 : {ownershipParticle1, ownedOrHalo}) {
            // skip inapplicable cases
            if (ownershipCell1 == owned and (ownershipParticle2 == halo or ownershipParticle3 == halo)) {
              continue;
            }
            if (ownershipCell1 == halo and (ownershipParticle2 == owned or ownershipParticle3 == owned)) {
              continue;
            }

            mdLib::AxilrodTellerMutoFunctor<Molecule> ATMFunctor(CellFunctorTest3B<CellFunctorType>::cutoff);
            ATMFunctor.setParticleProperties(CellFunctorTest3B<CellFunctorType>::nu);

            ATMFunctor.initTraversal();

            CellFunctorType cellFunctor(&ATMFunctor, CellFunctorTest3B<CellFunctorType>::cutoff);
            cellFunctor.setSortingThreshold(sortingThreshold);

            const auto &[forceParticle1, forceParticle2, forceParticle3] = doSingleCellInteraction<CellFunctorType>(
                cellFunctor, ownershipParticle1, ownershipParticle2, ownershipParticle3, ownershipCell1,
                cellFunctor.getDataLayout(), ATMFunctor);

            ATMFunctor.endTraversal(cellFunctor.getNewton3());

            // EXPECTATIONS
            if (ownershipCell1 == halo) {
              EXPECT_EQ(forceParticle1, 0) << createErrorString(ownershipCell1, ownershipParticle1, ownershipParticle2,
                                                                ownershipParticle3, 1, true);
              EXPECT_EQ(forceParticle2, 0) << createErrorString(ownershipCell1, ownershipParticle1, ownershipParticle2,
                                                                ownershipParticle3, 2, true);
              EXPECT_EQ(forceParticle3, 0) << createErrorString(ownershipCell1, ownershipParticle1, ownershipParticle2,
                                                                ownershipParticle3, 3, true);
            } else {
              // in all other cases we expect force on all particles in Cell1 as long as they are owned.
              if (ownershipParticle1 == owned) {
                EXPECT_GT(forceParticle1, 0) << createErrorString(ownershipCell1, ownershipParticle1,
                                                                  ownershipParticle2, ownershipParticle3, 1, false);
              }
              // if bidirectional or newton3=true we expect also force on particle 2
              if (ownershipParticle2 == owned) {
                EXPECT_GT(forceParticle2, 0) << createErrorString(ownershipCell1, ownershipParticle1,
                                                                  ownershipParticle2, ownershipParticle3, 2, false);
              }
              // if bidirectional or newton3=true we expect also force on particle 3
              if (ownershipParticle3 == owned) {
                EXPECT_GT(forceParticle3, 0) << createErrorString(ownershipCell1, ownershipParticle1,
                                                                  ownershipParticle2, ownershipParticle3, 3, false);
              }
            }
          }
        }
      }
    }
  }
}

/**
 * Tests the CellFunctor3B for a cell pair. Tests if force calculation is skipped if the base cell is a pure halo cell.
 * Checks if force calculations are done for a cell that can contain owned particles. Tests both cases where 2 particles
 * are in cell 1 and 1 particle in cell 2 as well as the other way around.
 */
TYPED_TEST_P(CellFunctorTest3B, testOwnedAndHaloCellInteractionPair) {
  using CellFunctorType = TypeParam;

  // shorthands for readability
  constexpr autopas::OwnershipState owned = autopas::OwnershipState::owned;
  constexpr autopas::OwnershipState halo = autopas::OwnershipState::halo;
  constexpr autopas::OwnershipState ownedOrHalo = autopas::OwnershipState::owned | autopas::OwnershipState::halo;

  // helper function to create a summarizing error string
  const auto createErrorString = [&](const auto ownershipCell1, const auto ownershipCell2, const auto ownershipP1,
                                     const auto ownershipP2, const bool particle2inCell1, const auto ownershipP3,
                                     const int numParticle, const bool actualForceZero) {
    std::stringstream errorStr;
    errorStr << "Particle " << numParticle << (actualForceZero ? " does not experience" : " experiences") << " force."
             << "\nCell1: " << ownershipCell1 << "\nCell2: " << ownershipCell2
             << "\nParticle 1 in Cell1: " << ownershipP1 << "\nParticle 2 in " << (particle2inCell1 ? "Cell1" : "Cell2")
             << ": " << ownershipP2 << "\nParticle 3 in Cell2: " << ownershipP3;
    return errorStr.str();
  };

  // Test with and without sorting
  for (const auto sortingThreshold : {0, 100}) {
    // Test all reasonable combinations of owned / halo particles and cells
    for (const auto ownershipParticle1 : {owned, halo}) {
      for (const auto ownershipParticle2 : {owned, halo}) {
        for (const auto ownershipParticle3 : {owned, halo}) {
          for (const auto ownershipCell1 : {ownershipParticle1, ownedOrHalo}) {
            for (const auto ownershipCell2 : {ownershipParticle3, ownedOrHalo}) {
              for (const auto particle2InCell1 : {true, false}) {
                // skip inapplicable cases
                if (ownershipCell1 == owned and particle2InCell1 and ownershipParticle2 == halo) {
                  continue;
                }
                if (ownershipCell1 == halo and particle2InCell1 and ownershipParticle2 == owned) {
                  continue;
                }
                if (ownershipCell2 == owned and not particle2InCell1 and ownershipParticle2 == halo) {
                  continue;
                }
                if (ownershipCell2 == halo and not particle2InCell1 and ownershipParticle2 == owned) {
                  continue;
                }
                mdLib::AxilrodTellerMutoFunctor<Molecule> ATMFunctor(CellFunctorTest3B<CellFunctorType>::cutoff);
                ATMFunctor.setParticleProperties(CellFunctorTest3B<CellFunctorType>::nu);

                ATMFunctor.initTraversal();

                CellFunctorType cellFunctor(&ATMFunctor, CellFunctorTest3B<CellFunctorType>::cutoff);
                cellFunctor.setSortingThreshold(sortingThreshold);

                const auto &[forceParticle1, forceParticle2, forceParticle3] = doPairCellInteraction<CellFunctorType>(
                    cellFunctor, ownershipParticle1, ownershipParticle2, ownershipParticle3, ownershipCell1,
                    ownershipCell2, cellFunctor.getDataLayout(), ATMFunctor, particle2InCell1);

                ATMFunctor.endTraversal(cellFunctor.getNewton3());

                // EXPECTATIONS
                if (ownershipCell1 == halo and ownershipCell2 == halo) {
                  EXPECT_EQ(forceParticle1, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                           particle2InCell1, ownershipParticle3, 1, true);
                  EXPECT_EQ(forceParticle2, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                           particle2InCell1, ownershipParticle3, 2, true);
                  EXPECT_EQ(forceParticle3, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                           particle2InCell1, ownershipParticle3, 3, true);
                } else if (ownershipCell1 == halo and (not cellFunctor.getNewton3()) and
                           (not cellFunctor.getBidirectional())) {
                  EXPECT_EQ(forceParticle1, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                           particle2InCell1, ownershipParticle3, 1, true);
                  if (particle2InCell1) {
                    EXPECT_EQ(forceParticle2, 0)
                        << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                             particle2InCell1, ownershipParticle3, 2, true);
                  }
                } else {
                  // in all other cases we expect force on particle A in Cell1 as long as it is owned.
                  if (ownershipParticle1 == owned) {
                    EXPECT_GT(forceParticle1, 0)
                        << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                             particle2InCell1, ownershipParticle3, 1, false);
                  }
                  // if bidirectional or newton3=true we expect also force on particle 2
                  if (ownershipParticle2 == owned and
                      (particle2InCell1 or (cellFunctor.getNewton3() or cellFunctor.getBidirectional()))) {
                    EXPECT_GT(forceParticle2, 0)
                        << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                             particle2InCell1, ownershipParticle3, 2, false);
                  }
                  // if bidirectional or newton3=true we expect also force on particle 3
                  if (ownershipParticle3 == owned and (cellFunctor.getNewton3() or cellFunctor.getBidirectional())) {
                    EXPECT_GT(forceParticle3, 0)
                        << createErrorString(ownershipCell1, ownershipCell2, ownershipParticle1, ownershipParticle2,
                                             particle2InCell1, ownershipParticle3, 3, false);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

/**
 * Tests the CellFunctor3B for a cell triplet. Tests if force calculation is skipped if the base cell is a pure halo
 * cell. Checks if force calculations are done for a cell that can contain owned particles.
 */
TYPED_TEST_P(CellFunctorTest3B, testOwnedAndHaloCellInteractionTriple) {
  using CellFunctorType = TypeParam;

  // shorthands for readability
  constexpr autopas::OwnershipState owned = autopas::OwnershipState::owned;
  constexpr autopas::OwnershipState halo = autopas::OwnershipState::halo;
  constexpr autopas::OwnershipState ownedOrHalo = autopas::OwnershipState::owned | autopas::OwnershipState::halo;

  // helper function to create a summarizing error string
  const auto createErrorString = [&](const auto ownershipCell1, const auto ownershipCell2, const auto ownershipCell3,
                                     const auto ownershipP1, const auto ownershipP2, const auto ownershipP3,
                                     const int numParticle, const bool actualForceZero) {
    std::stringstream errorStr;
    errorStr << "Particle " << numParticle << (actualForceZero ? " does not experience" : " experiences") << " force."
             << "\nCell1: " << ownershipCell1 << "\nCell2: " << ownershipCell2 << "\nCell3: " << ownershipCell3
             << "\nParticle 1 in Cell1: " << ownershipP1 << "\nParticle 2 in Cell2: " << ownershipP2
             << "\nParticle 3 in Cell3: " << ownershipP3;
    return errorStr.str();
  };

  // Test with and without sorting
  for (const auto sortingThreshold : {0, 100}) {
    // Test all reasonable combinations of owned / halo particles and cells
    for (const auto ownershipParticle1 : {owned, halo}) {
      for (const auto ownershipParticle2 : {owned, halo}) {
        for (const auto ownershipParticle3 : {owned, halo}) {
          for (const auto ownershipCell1 : {ownershipParticle1, ownedOrHalo}) {
            for (const auto ownershipCell2 : {ownershipParticle2, ownedOrHalo}) {
              for (const auto ownershipCell3 : {ownershipParticle3, ownedOrHalo}) {
                mdLib::AxilrodTellerMutoFunctor<Molecule> ATMFunctor(CellFunctorTest3B<CellFunctorType>::cutoff);
                ATMFunctor.setParticleProperties(CellFunctorTest3B<CellFunctorType>::nu);

                ATMFunctor.initTraversal();

                CellFunctorType cellFunctor(&ATMFunctor, CellFunctorTest3B<CellFunctorType>::cutoff);
                cellFunctor.setSortingThreshold(sortingThreshold);

                const auto &[forceParticle1, forceParticle2, forceParticle3] = doTripleCellInteraction<CellFunctorType>(
                    cellFunctor, ownershipParticle1, ownershipParticle2, ownershipParticle3, ownershipCell1,
                    ownershipCell2, ownershipCell3, cellFunctor.getDataLayout(), ATMFunctor);

                ATMFunctor.endTraversal(cellFunctor.getNewton3());

                // EXPECTATIONS
                if (ownershipCell1 == halo and ownershipCell2 == halo and ownershipCell3 == halo) {
                  EXPECT_EQ(forceParticle1, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipCell3, ownershipParticle1,
                                           ownershipParticle2, ownershipParticle3, 1, true);
                  EXPECT_EQ(forceParticle2, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipCell3, ownershipParticle1,
                                           ownershipParticle2, ownershipParticle3, 2, true);
                  EXPECT_EQ(forceParticle3, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipCell3, ownershipParticle1,
                                           ownershipParticle2, ownershipParticle3, 3, true);
                } else if (ownershipCell1 == halo and (not cellFunctor.getNewton3()) and
                           (not cellFunctor.getBidirectional())) {
                  EXPECT_EQ(forceParticle1, 0)
                      << createErrorString(ownershipCell1, ownershipCell2, ownershipCell3, ownershipParticle1,
                                           ownershipParticle2, ownershipParticle3, 1, true);
                } else {
                  // in all other cases we expect force on particle A in Cell1 as long as it is owned.
                  if (ownershipParticle1 == owned) {
                    EXPECT_GT(forceParticle1, 0)
                        << createErrorString(ownershipCell1, ownershipCell2, ownershipCell3, ownershipParticle1,
                                             ownershipParticle2, ownershipParticle3, 1, false);
                  }
                  // if bidirectional or newton3=true we expect also force on particle 2
                  if (ownershipParticle2 == owned and (cellFunctor.getNewton3() or cellFunctor.getBidirectional())) {
                    EXPECT_GT(forceParticle2, 0)
                        << createErrorString(ownershipCell1, ownershipCell2, ownershipCell3, ownershipParticle1,
                                             ownershipParticle2, ownershipParticle3, 2, false);
                  }
                  // if bidirectional or newton3=true we expect also force on particle 3
                  if (ownershipParticle3 == owned and (cellFunctor.getNewton3() or cellFunctor.getBidirectional())) {
                    EXPECT_GT(forceParticle3, 0)
                        << createErrorString(ownershipCell1, ownershipCell2, ownershipCell3, ownershipParticle1,
                                             ownershipParticle2, ownershipParticle3, 3, false);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(CellFunctorTest3B, testOwnedAndHaloCellInteractionSingle,
                            testOwnedAndHaloCellInteractionPair, testOwnedAndHaloCellInteractionTriple);
INSTANTIATE_TYPED_TEST_SUITE_P(TypedTest, CellFunctorTest3B, CellFTestingTypes);