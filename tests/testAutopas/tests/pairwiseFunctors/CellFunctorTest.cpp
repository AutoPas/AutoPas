/**
 * @file CellFunctorTest.cpp
 * @author D. Martin
 * @date 29.08.23
 */

#include "CellFunctorTest.h"

// Type aliases via inheritance for more readable test names (using declarations do not work for this)
struct CellFunctor_AoS_NoN3_NoBi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::aos, false, false> {
  CellFunctor_AoS_NoN3_NoBi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor(f, sortingCutoff) {}
};
struct CellFunctor_AoS_NoN3_Bi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::aos, false, true> {
  CellFunctor_AoS_NoN3_Bi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff) : CellFunctor(f, sortingCutoff) {}
};
struct CellFunctor_AoS_N3_NoBi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::aos, true, false> {
  CellFunctor_AoS_N3_NoBi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff) : CellFunctor(f, sortingCutoff) {}
};
struct CellFunctor_AoS_N3_Bi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::aos, true, true> {
  CellFunctor_AoS_N3_Bi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff) : CellFunctor(f, sortingCutoff) {}
};
struct CellFunctor_SoA_NoN3_NoBi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::soa, false, false> {
  CellFunctor_SoA_NoN3_NoBi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff)
      : CellFunctor(f, sortingCutoff) {}
};
struct CellFunctor_SoA_NoN3_Bi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::soa, false, true> {
  CellFunctor_SoA_NoN3_Bi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff) : CellFunctor(f, sortingCutoff) {}
};
struct CellFunctor_SoA_N3_NoBi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::soa, true, false> {
  CellFunctor_SoA_N3_NoBi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff) : CellFunctor(f, sortingCutoff) {}
};
struct CellFunctor_SoA_N3_Bi
    : public autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                            autopas::DataLayoutOption::soa, true, true> {
  CellFunctor_SoA_N3_Bi(mdLib::LJFunctor<Molecule> *f, const double sortingCutoff) : CellFunctor(f, sortingCutoff) {}
};

/**
 * All relevant CellFunctor configurations which are used to instantiate the typed test cases
 * testOwnedAndHaloCellInteractionPair and testOwnedAndHaloCellInteractionSingle
 *
 */
// clang-format off
using CellFTestingTypes = ::testing::Types<CellFunctor_AoS_NoN3_NoBi,
                                           CellFunctor_AoS_NoN3_Bi,
                                           CellFunctor_AoS_N3_NoBi,
                                           CellFunctor_AoS_N3_Bi,
                                           CellFunctor_SoA_NoN3_NoBi,
                                           CellFunctor_SoA_NoN3_Bi,
                                           CellFunctor_SoA_N3_NoBi,
                                           CellFunctor_SoA_N3_Bi>;
// clang-format on

/**
 * Helper function for readability
 * @param ownershipState
 * @return True if the passed ownership state contains the bit for "owned".
 */
bool containsOwned(autopas::OwnershipState ownershipState) {
  return static_cast<bool>(ownershipState & autopas::OwnershipState::owned);
}

/**
 * Helper for testOwnedAndHaloCellInteractionPair andtestOwnedAndHaloCellInteractionSingle that creates particles and
 * cells with given OwnershipStates, executes CellFunctor::processCellPair or CellFunctor::processCell and returns the
 * calculated forces on particles.
 *
 * @tparam T The type of CellFunctor
 * @param cellFunctor
 * @param ownershipParticle1 OwnershipState for particle 1 (for testOwnedAndHaloCellInteractionPair this particle is in
 * cell 1)
 * @param ownershipParticle2 OwnershipState for particle 2 (for testOwnedAndHaloCellInteractionPair this particle is in
 * cell 2)
 * @param ownershipCell1 OwnershipState for cell 1 (for testOwnedAndHaloCellInteractionSingle both particles are in this
 * cell)
 * @param ownershipCell2 OwnershipState for cell 2 (only used in testOwnedAndHaloCellInteractionPair)
 * @param dataLayout AoS or SoA
 * @param ljFunctor The functor to evaluate forces
 * @param singleCell boolean if called from testOwnedAndHaloCellInteractionPair (false) or from
 * testOwnedAndHaloCellInteractionSingle (true)
 * @return The force in X direction of the first particle in the first and second cell. In the single cell case
 * the same force is returned twice. std::tuple<f0[0], f1[0]>
 */
template <class T>
std::tuple<double, double> ownedHaloInteractionHelper(T &cellFunctor, const autopas::OwnershipState ownershipParticle1,
                                                      const autopas::OwnershipState ownershipParticle2,
                                                      const autopas::OwnershipState ownershipCell1,
                                                      const autopas::OwnershipState ownershipCell2,
                                                      const autopas::DataLayoutOption dataLayout,
                                                      mdLib::LJFunctor<Molecule> &ljFunctor, bool singleCell) {
  // create two particles
  Molecule p1({0.6, 0.5, 0.5}, {0., 0., 0.}, 0);
  p1.setOwnershipState(ownershipParticle1);

  Molecule p2({1.4, 0.5, 0.5}, {0., 0., 0.}, 1);
  p2.setOwnershipState(ownershipParticle2);

  // create cells
  autopas::FullParticleCell<Molecule> cell1({1., 1., 1.});
  autopas::FullParticleCell<Molecule> cell2({1., 1., 1.});

  cell1.setPossibleParticleOwnerships(ownershipCell1);
  cell1.addParticle(p1);

  if (not singleCell) {
    cell2.setPossibleParticleOwnerships(ownershipCell2);
    cell2.addParticle(p2);
  } else {
    cell1.addParticle(p2);
  }

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ljFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
    if (not singleCell) {
      ljFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
    }
  }

  if (not singleCell) {
    cellFunctor.processCellPair(cell1, cell2);
  } else {
    cellFunctor.processCell(cell1);
  }

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ljFunctor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    if (not singleCell) {
      ljFunctor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    }
  }

  if (not singleCell) {
    return std::tuple<double, double>{std::abs(cell1[0].getF()[0]), std::abs(cell2[0].getF()[0])};
  } else {
    return std::tuple<double, double>{std::abs(cell1[0].getF()[0]), std::abs(cell1[1].getF()[0])};
  }
}

TYPED_TEST_SUITE_P(CellFunctorTest);

/**
 *
 * Tests if the cell functor skips the cell interactions:
 *  - halo <-> halo
 *  - halo  -> any
 *
 * All other interaction combinations should be applied.
 */
TYPED_TEST_P(CellFunctorTest, testOwnedAndHaloCellInteractionPair) {
  using CellFunctorType = TypeParam;

  constexpr double cutoff = 1.;
  constexpr double sigma = 1.;
  constexpr double epsilon = 1.;

  // shorthands for readability
  constexpr autopas::OwnershipState owned = autopas::OwnershipState::owned;
  constexpr autopas::OwnershipState halo = autopas::OwnershipState::halo;
  constexpr autopas::OwnershipState ownedOrHalo = autopas::OwnershipState::owned | autopas::OwnershipState::halo;

  // Test all reasonable combinations of owned / halo particles and cells
  for (const auto ownershipParticleA : {owned, halo}) {
    for (const auto ownershipParticleB : {owned, halo}) {
      for (const auto ownershipCellA : {ownershipParticleA, ownedOrHalo}) {
        for (const auto ownerShipStateCellB : {ownershipParticleB, ownedOrHalo}) {
          mdLib::LJFunctor<Molecule> ljFunctor(cutoff);
          ljFunctor.setParticleProperties(sigma, epsilon);

          ljFunctor.initTraversal();

          CellFunctorType cellFunctor(&ljFunctor, cutoff);

          const auto &[forceParticleA, forceParticleB] = ownedHaloInteractionHelper<CellFunctorType>(
              cellFunctor, ownershipParticleA, ownershipParticleB, ownershipCellA, ownerShipStateCellB,
              cellFunctor.getDataLayout(), ljFunctor, false);

          ljFunctor.endTraversal(cellFunctor.getNewton3());

          // EXPECTATIONS
          if (ownershipCellA == halo and ownerShipStateCellB == halo) {
            EXPECT_EQ(forceParticleA, 0)
                << "Particle in Cell1 is experiencing force with both cells containing only halo particles."
                   "\nCell1: "
                << ownershipCellA << "\nCell2: " << ownerShipStateCellB << "\nParticle in Cell1: " << ownershipParticleA
                << "\nParticle in Cell2: " << ownershipParticleB;
            EXPECT_EQ(forceParticleB, 0)
                << "Particle in Cell2 is experiencing force with with both cells containing only halo particles."
                   "\nCell1: "
                << ownershipCellA << "\nCell2: " << ownerShipStateCellB << "\nParticle in Cell1: " << ownershipParticleA
                << "\nParticle in Cell2: " << ownershipParticleB;
          } else if (ownershipCellA == halo and (not cellFunctor.getNewton3()) and
                     (not cellFunctor.getBidirectional())) {
            // if cell1 is halo, and NoN3 and no bidirectional we can skip the interaction
            EXPECT_EQ(forceParticleA, 0) << "Particle in Cell1 is experiencing force with "
                                            "OwnershipState halo, no newton3 and bidirectional off."
                                            "\nCell1: "
                                         << ownershipCellA << "\nCell2: " << ownerShipStateCellB
                                         << "\nParticle in Cell1: " << ownershipParticleA
                                         << "\nParticle in Cell2: " << ownershipParticleB;
          } else {
            // in all other cases we expect force on particle in Cell1
            EXPECT_GT(forceParticleA, 0) << "Particle in Cell1 does not experience force."
                                            "\nCell1: "
                                         << ownershipCellA << "\nCell2: " << ownerShipStateCellB
                                         << "\nParticle in Cell1: " << ownershipParticleA
                                         << "\nParticle in Cell2: " << ownershipParticleB;

            // if bidirectional or newton3=true we expect also force on particle in Cell2
            if ((cellFunctor.getBidirectional() or cellFunctor.getNewton3()) and containsOwned(ownershipParticleB)) {
              EXPECT_GT(forceParticleB, 0)
                  << "Particle in Cell2 does not experience force."
                     "\nCell1: "
                  << ownershipCellA << "\nCell2: " << ownerShipStateCellB
                  << "\nParticle in Cell1: " << ownershipParticleA << "\nParticle in Cell2: " << ownershipParticleB;
            }
          }
        }
      }
    }
  }
}

/**
 * Tests if force calculation is skipped in the CellFunctor with a pure halo cell. Checks if force calculations are
 * done for a cell that can contain owned particles.
 *
 */
TYPED_TEST_P(CellFunctorTest, testOwnedAndHaloCellInteractionSingle) {
  using CellFunctorType = TypeParam;

  constexpr double cutoff = 1.;
  constexpr double sigma = 1.;
  constexpr double epsilon = 1.;

  // shorthands for readability
  constexpr autopas::OwnershipState owned = autopas::OwnershipState::owned;
  constexpr autopas::OwnershipState halo = autopas::OwnershipState::halo;
  constexpr autopas::OwnershipState ownedOrHalo = autopas::OwnershipState::owned | autopas::OwnershipState::halo;

  // Test all reasonable combinations of owned / halo particles and cells
  for (const auto ownershipParticleA : {owned, halo}) {
    for (const auto ownershipParticleB : {owned, halo}) {
      for (const auto ownershipCellA : {ownershipParticleA, ownedOrHalo}) {
        // skip inapplicable cases
        if (ownershipCellA == owned and ownershipParticleB == halo) {
          continue;
        }
        if (ownershipCellA == halo and ownershipParticleB == owned) {
          continue;
        }

        mdLib::LJFunctor<Molecule> ljFunctor(cutoff);
        ljFunctor.setParticleProperties(sigma, epsilon);

        ljFunctor.initTraversal();

        CellFunctorType cellFunctor(&ljFunctor, cutoff);

        const auto &[forceParticleA, forceParticleB] = ownedHaloInteractionHelper<CellFunctorType>(
            cellFunctor, ownershipParticleA, ownershipParticleB, ownershipCellA, ownershipCellA,
            cellFunctor.getDataLayout(), ljFunctor, true);

        ljFunctor.endTraversal(cellFunctor.getNewton3());

        // EXPECTATIONS
        if (ownershipCellA == halo) {
          EXPECT_EQ(forceParticleA, 0) << "Particle 1 is experiencing force."
                                          "\nCell1: "
                                       << ownershipCellA << "\nParticle 1 in Cell1: " << ownershipParticleA
                                       << "\nParticle 2 in Cell1: " << ownershipParticleB;
          EXPECT_EQ(forceParticleB, 0) << "Particle 2 is experiencing force."
                                          "\nCell1: "
                                       << ownershipCellA << "\nParticle 1 in Cell1: " << ownershipParticleA
                                       << "\nParticle 2 in Cell1: " << ownershipParticleB;
        } else {
          // in all other cases we expect force on particle A in Cell1 as long as it is owned.
          if (containsOwned(ownershipParticleA)) {
            EXPECT_GT(forceParticleA, 0) << "Particle 1 in Cell1 does not experience force."
                                            "\nCell1: "
                                         << ownershipCellA << "\nParticle 1 in Cell1: " << ownershipParticleA
                                         << "\nParticle 2 in Cell1: " << ownershipParticleB;
          }
          // if bidirectional or newton3=true we expect also force on particle in Cell2
          if ((cellFunctor.getBidirectional() or cellFunctor.getNewton3()) and containsOwned(ownershipParticleB)) {
            EXPECT_GT(forceParticleB, 0) << "Particle 2 does not experience force."
                                            "\nCell1: "
                                         << ownershipCellA << "\nParticle 1 in Cell1: " << ownershipParticleA
                                         << "\nParticle 2 in Cell1: " << ownershipParticleB;
          }
        }
      }
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(CellFunctorTest, testOwnedAndHaloCellInteractionPair,
                            testOwnedAndHaloCellInteractionSingle);
INSTANTIATE_TYPED_TEST_SUITE_P(TypedTest, CellFunctorTest, CellFTestingTypes);