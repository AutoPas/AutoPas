/**
 * @file CellFunctorTest.cpp
 * @author D. Martin
 * @date 29.08.23
 */

#include "CellFunctorTest.h"

/**
 * All relevant CellFunctor configurations which are used to instantiate the typed test cases
 * testOwnedAndHaloCellInteractionPair and testOwnedAndHaloCellInteractionSingle
 *
 */
using CellFTestingTypes = ::testing::Types<
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, false, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, false, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, false, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, false, true>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, true, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, true, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, true, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, true, true>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, false, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, false, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, false, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, false, true>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, true, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, true, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, true, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, true, true>>>;

/**
 * Helper for testOwnedAndHaloCellInteractionPair andtestOwnedAndHaloCellInteractionSingle that creates particles and
 * cells with given OwnershipStates, executes CellFunctor::processCellPair or CellFunctor::processCell and returnes the
 * calculated forces on particles.
 *
 * @tparam T The type of CellFunctor
 * @param cellFunctor
 * @param osp1 OwnershipState for particle 1 (for testOwnedAndHaloCellInteractionPair this partile is in cell 1)
 * @param osp2 OwnershipState for particle 2 (for testOwnedAndHaloCellInteractionPair this partile is in cell 2)
 * @param osc1 OwnershipState for cell 1 (for testOwnedAndHaloCellInteractionSingle both particles are in this cell)
 * @param osc2 OwnershipState for cell 2 (only used in testOwnedAndHaloCellInteractionPair)
 * @param dataLayout AoS or SoA
 * @param ljFunctor The functor to evaluate forces
 * @param singleCell boolean if called from testOwnedAndHaloCellInteractionPair (false) or from
 * testOwnedAndHaloCellInteractionSingle (true)
 * @return std::tuple<double, double>
 */
template <class T>
std::tuple<double, double> ownedHaloInteractionHelper(T &cellFunctor, const autopas::OwnershipState osp1,
                                                      const autopas::OwnershipState osp2,
                                                      const autopas::OwnershipState osc1,
                                                      const autopas::OwnershipState osc2,
                                                      const autopas::DataLayoutOption dataLayout,
                                                      mdLib::LJFunctor<Molecule> &ljFunctor, bool singleCell) {
  // create two particles
  Molecule p1({0.6, 0.5, 0.5}, {0., 0., 0.}, 0);
  p1.setOwnershipState(osp1);

  Molecule p2({1.4, 0.5, 0.5}, {0., 0., 0.}, 1);
  p2.setOwnershipState(osp2);

  // create cells
  autopas::FullParticleCell<Molecule> cell1({1., 1., 1.});
  autopas::FullParticleCell<Molecule> cell2({1., 1., 1.});

  cell1.setPossibleParticleOwnerships(osc1);
  cell1.addParticle(p1);

  if (not singleCell) {
    cell2.setPossibleParticleOwnerships(osc2);
    cell2.addParticle(p2);
  } else {
    cell1.addParticle(p2);
  }

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ljFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0);
    if (not singleCell) {
      ljFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0);
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

TYPED_TEST_CASE_P(CellFunctorTest);

/**
 * Tests if pure halo-halo cell interactions or interactions between a halo cell and any other cell with newton3==false
 * and bidirection==false are skipped in the CellFunctor (no forces are calculated). Tests if all other interactions
 * result in force calculations.
 *
 */
TYPED_TEST_P(CellFunctorTest, testOwnedAndHaloCellInteractionPair) {
  const double cutoff = 1.;
  const double sigma = 1.;
  const double epsilon = 1.;

  using CellFT = typename TypeParam::first_type;
  using CFWT = typename TypeParam::second_type;

  for (auto ownerShipStateA : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    for (auto ownerShipStateB : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
      for (auto ownerShipStateCellA :
           {ownerShipStateA, (autopas::OwnershipState::owned | autopas::OwnershipState::halo)}) {
        for (auto ownerShipStateCellB :
             {ownerShipStateB, (autopas::OwnershipState::owned | autopas::OwnershipState::halo)}) {
          mdLib::LJFunctor<Molecule> ljFunctor(cutoff);
          ljFunctor.setParticleProperties(sigma, epsilon);

          ljFunctor.initTraversal();

          CellFT cellFunctor(&ljFunctor, cutoff);
          CFWT cfw;

          std::tuple<double, double> result =
              ownedHaloInteractionHelper<CellFT>(cellFunctor, ownerShipStateA, ownerShipStateB, ownerShipStateCellA,
                                                 ownerShipStateCellB, cfw.dataLayoutV, ljFunctor, false);

          ljFunctor.endTraversal(cfw.useNewton3V);

          if (ownerShipStateCellA == autopas::OwnershipState::halo and
              ownerShipStateCellB == autopas::OwnershipState::halo) {
            EXPECT_NEAR(std::get<0>(result), 0, 1e-13)
                << "Particle in Cell1 is experiencing force with both cells containing only halo particles. "
                   "OwnershipState of Cell1: "
                << ownerShipStateCellA << ", OwnershipState of Cell2: " << ownerShipStateCellB
                << ", OwnershipState of Particle in Cell1: " << ownerShipStateA
                << ", OwnershipState of Particle in Cell2: " << ownerShipStateB;
            EXPECT_NEAR(std::get<1>(result), 0, 1e-13)
                << "Particle in Cell2 is experiencing force with with both cells containing only halo particles. "
                   "OwnershipState of Cell1: "
                << ownerShipStateCellA << ", OwnershipState of Cell2: " << ownerShipStateCellB
                << ", OwnershipState of Particle in Cell1: " << ownerShipStateA
                << ", OwnershipState of Particle in Cell2: " << ownerShipStateB;
          } else if (ownerShipStateCellA == autopas::OwnershipState::halo and (not cfw.useNewton3V) and
                     (not cfw.bidirectionalV)) {
            // if cell1 is halo, and NoN3 and no bidirectional we can skip the interaction
            EXPECT_NEAR(std::get<0>(result), 0, 1e-13)
                << "Particle in Cell1 is experiencing force with "
                   "OwnershipState halo, no newton3 and bidirectional off. "
                   "OwnershipState of Cell1: "
                << ownerShipStateCellA << ", OwnershipState of Cell2: " << ownerShipStateCellB
                << ", OwnershipState of Particle in Cell1: " << ownerShipStateA
                << ", OwnershipState of Particle in Cell2: " << ownerShipStateB;
          } else {
            // in all other cases we expect force on particle in Cell1
            EXPECT_TRUE(std::get<0>(result) > 0)
                << "Particle in Cell1 does not experience force. OwnershipState of Cell1: " << ownerShipStateCellA
                << ", OwnershipState of Cell2: " << ownerShipStateCellB
                << ", OwnershipState of Particle in Cell1: " << ownerShipStateA
                << ", OwnershipState of Particle in Cell2: " << ownerShipStateB;

            // if bidirectional or newton3=true we expect also force on particle in Cell2
            if (cfw.bidirectionalV or cfw.useNewton3V) {
              EXPECT_TRUE(std::get<1>(result) > 0)
                  << "Particle in Cell2 does not experience force. OwnershipState of Cell1: " << ownerShipStateCellA
                  << ", OwnershipState of Cell2: " << ownerShipStateCellB
                  << ", OwnershipState of Particle in Cell1: " << ownerShipStateA
                  << ", OwnershipState of Particle in Cell2: " << ownerShipStateB;
            }
          }
        }
      }
    }
  }
}

/**
 * Tests if force calculation is skipped in the CellFunctor with a pure halo cell. Checks if force calculations are done
 * for a cell that can contain owned particles.
 *
 */
TYPED_TEST_P(CellFunctorTest, testOwnedAndHaloCellInteractionSingle) {
  const double cutoff = 1.;
  const double sigma = 1.;
  const double epsilon = 1.;

  using CellFT = typename TypeParam::first_type;
  using CFWT = typename TypeParam::second_type;

  for (auto ownerShipStateA : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    for (auto ownerShipStateB : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
      for (auto ownerShipStateCellA :
           {ownerShipStateA, (autopas::OwnershipState::owned | autopas::OwnershipState::halo)}) {
        // skip unapplicable cases
        if (ownerShipStateCellA == autopas::OwnershipState::owned and
            ownerShipStateB == autopas::OwnershipState::halo) {
          continue;
        }
        if (ownerShipStateCellA == autopas::OwnershipState::halo and
            ownerShipStateB == autopas::OwnershipState::owned) {
          continue;
        }

        mdLib::LJFunctor<Molecule> ljFunctor(cutoff);
        ljFunctor.setParticleProperties(sigma, epsilon);

        ljFunctor.initTraversal();

        CellFT cellFunctor(&ljFunctor, cutoff);
        CFWT cfw;

        std::tuple<double, double> result =
            ownedHaloInteractionHelper<CellFT>(cellFunctor, ownerShipStateA, ownerShipStateB, ownerShipStateCellA,
                                               ownerShipStateCellA, cfw.dataLayoutV, ljFunctor, true);

        ljFunctor.endTraversal(cfw.useNewton3V);

        if (ownerShipStateCellA == autopas::OwnershipState::halo) {
          EXPECT_NEAR(std::get<0>(result), 0, 1e-13)
              << "Particle 1 is experiencing force. OwnershipState of Cell1: " << ownerShipStateCellA
              << ", OwnershipState of Particle 1 in Cell1: " << ownerShipStateA
              << ", OwnershipState of Particle 2 in Cell1: " << ownerShipStateB;
          EXPECT_NEAR(std::get<1>(result), 0, 1e-13)
              << "Particle 2 is experiencing force. OwnershipState of Cell1: " << ownerShipStateCellA
              << ", OwnershipState of Particle 1 in Cell1: " << ownerShipStateA
              << ", OwnershipState of Particle 2 in Cell1: " << ownerShipStateB;
        } else {
          // in all other cases we expect force on particle in Cell1
          EXPECT_TRUE(std::get<0>(result) > 0)
              << "Particle 1 in Cell1 does not experience force. OwnershipState of Cell1: " << ownerShipStateCellA
              << ", OwnershipState of Particle 1 in Cell1: " << ownerShipStateA
              << ", OwnershipState of Particle 2 in Cell1: " << ownerShipStateB;

          // if bidirectional or newton3=true we expect also force on particle in Cell2
          if (cfw.bidirectionalV or cfw.useNewton3V) {
            EXPECT_TRUE(std::get<1>(result) > 0)
                << "Particle 2 does not experience force. OwnershipState of Cell1: " << ownerShipStateCellA
                << ", OwnershipState of Particle 1 in Cell1: " << ownerShipStateA
                << ", OwnershipState of Particle 2 in Cell1: " << ownerShipStateB;
          }
        }
      }
    }
  }
}

REGISTER_TYPED_TEST_CASE_P(CellFunctorTest, testOwnedAndHaloCellInteractionPair, testOwnedAndHaloCellInteractionSingle);
INSTANTIATE_TYPED_TEST_CASE_P(TypedTest, CellFunctorTest, CellFTestingTypes);