/**
 * @file CellOwnershipTest.cpp
 * @author D. Martin
 * @date 10.09.23
 */

#include "CellOwnershipTest.h"

/**
 * Note: For OctreeNodeWrapper and ClusterTowers (which also inherit from ParticleCell) the tests in this file are not
 * executed.
 * The OctreeNodeWrapper only manages internal nodes. Only the OctreeLeafNodes store the particles in an internal
 * storage.
 * The ClusterTowers internally use a FullParticleCell to store particles, therefore these tests are not executed for
 * ClusterTowers, since it is covered by the tests for FullParticleCell.
 */

/**
 * Checks the default OwnershipState of FullParticleCell and expected exception behavior if OwnershipState is reset.
 * Tries to add owned and halo particles to cells with different ownership states and checks the expected exception
 * behavior.
 */
TEST(CellOwnershipTest, testOwnershipFullParticleCell) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::FullParticleCell<Molecule> cell({1., 1., 1.});

    // expect OwnershipState ownedOrHalo
    EXPECT_TRUE(cell.getPossibleParticleOwnerships() ==
                (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

    cell.setPossibleParticleOwnerships(ownershipState);

    // cell ownership can not be changed once set
    EXPECT_ANY_THROW(cell.setPossibleParticleOwnerships(ownershipState));

    Molecule pOwned({0.1, 0.1, 0.1}, {0., 0., 0.}, 0);
    pOwned.setOwnershipState(autopas::OwnershipState::owned);

    Molecule pHalo({0.2, 0.2, 0.2}, {0., 0., 0.}, 0);
    pHalo.setOwnershipState(autopas::OwnershipState::halo);

    Molecule pDummy({0.3, 0.3, 0.3}, {0., 0., 0.}, 0);
    pDummy.setOwnershipState(autopas::OwnershipState::dummy);

    Molecule pGood, pBad;

    if (ownershipState == autopas::OwnershipState::owned) {
      pGood = pOwned;
      pBad = pHalo;
    } else {
      pGood = pHalo;
      pBad = pOwned;
    }

    EXPECT_NO_THROW(cell.addParticle(pGood));

    // A dummy can be added always
    EXPECT_NO_THROW(cell.addParticle(pDummy));

    // A halo particle can not be added to a pure owned cell and vice versa
    EXPECT_ANY_THROW(cell.addParticle(pBad));
  }
}

/**
 * Checks the default OwnershipState of ReferenceParticleCell and expected exception behavior if OwnershipState is
 * reset.
 * Tries to add owned and halo particles to cells with different ownership states and checks the expected
 * exception behavior.
 */
TEST(CellOwnershipTest, testOwnershipReferenceParticleCell) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::ReferenceParticleCell<Molecule> cell({1., 1., 1.});

    // expect OwnershipState ownedOrHalo
    EXPECT_TRUE(cell.getPossibleParticleOwnerships() ==
                (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

    cell.setPossibleParticleOwnerships(ownershipState);

    // cell ownership can not be changed once set
    EXPECT_ANY_THROW(cell.setPossibleParticleOwnerships(ownershipState));

    Molecule pOwned({0.1, 0.1, 0.1}, {0., 0., 0.}, 0);
    pOwned.setOwnershipState(autopas::OwnershipState::owned);

    Molecule pHalo({0.2, 0.2, 0.2}, {0., 0., 0.}, 0);
    pHalo.setOwnershipState(autopas::OwnershipState::halo);

    Molecule pDummy({0.3, 0.3, 0.3}, {0., 0., 0.}, 0);
    pDummy.setOwnershipState(autopas::OwnershipState::dummy);

    Molecule pGood, pBad;

    if (ownershipState == autopas::OwnershipState::owned) {
      pGood = pOwned;
      pBad = pHalo;
    } else {
      pGood = pHalo;
      pBad = pOwned;
    }

    EXPECT_NO_THROW(cell.addParticleReference(&pGood));

    // A dummy can be added always
    EXPECT_NO_THROW(cell.addParticleReference(&pDummy));

    // A halo particle can not be added to a pure owned cell and vice versa
    EXPECT_ANY_THROW(cell.addParticleReference(&pBad));
  }
}

/**
 * Checks the default OwnershipState of OctreeLeafNode and expected exception behavior if OwnershipState is reset.
 * Tries to add owned and halo particles to cells with different ownership states and checks the expected exception
 * behavior.
 */
TEST(CellOwnershipTest, testOwnershipOctreeLeafNode) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::OctreeLeafNode<Molecule> cell({0, 0, 0}, {1., 1., 1.}, nullptr, 5, 1., 1.);

    // expect OwnershipState ownedOrHalo
    EXPECT_TRUE(cell.getPossibleParticleOwnerships() ==
                (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

    cell.setPossibleParticleOwnerships(ownershipState);

    EXPECT_TRUE(cell.getPossibleParticleOwnerships() == ownershipState);

    EXPECT_ANY_THROW(cell.setPossibleParticleOwnerships(ownershipState));

    Molecule pOwned({0.1, 0.1, 0.1}, {0., 0., 0.}, 0);
    pOwned.setOwnershipState(autopas::OwnershipState::owned);

    Molecule pHalo({0.2, 0.2, 0.2}, {0., 0., 0.}, 0);
    pHalo.setOwnershipState(autopas::OwnershipState::halo);

    Molecule pDummy({0.3, 0.3, 0.3}, {0., 0., 0.}, 0);
    pDummy.setOwnershipState(autopas::OwnershipState::dummy);

    Molecule pGood, pBad;

    if (ownershipState == autopas::OwnershipState::owned) {
      pGood = pOwned;
      pBad = pHalo;
    } else {
      pGood = pHalo;
      pBad = pOwned;
    }

    EXPECT_NO_THROW(cell.insert(pGood));

    // A dummy can be added always
    EXPECT_NO_THROW(cell.insert(pDummy));

    // A halo particle can not be added to a pure owned cell and vice versa
    EXPECT_ANY_THROW(cell.insert(pBad));
  }
}

/**
 * Tries to create a SortedCellView from a FullParticleCell and checks if the OwnershipState is transfered correctly
 * Note: Since a SortedCellView does not allow to add particles we do not check the exception behavior for adding
 * particles here.
 */
TEST(CellOwnershipTest, testOwnershipSortedCellView) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::FullParticleCell<Molecule> cell({1., 1., 1.});

    cell.setPossibleParticleOwnerships(ownershipState);

    auto scv = autopas::SortedCellView<Molecule, autopas::FullParticleCell<Molecule>>(cell, {0., 0., 0.});

    EXPECT_TRUE(scv.getPossibleParticleOwnerships() == ownershipState);

    EXPECT_ANY_THROW(scv.setPossibleParticleOwnerships(ownershipState));
  }
}