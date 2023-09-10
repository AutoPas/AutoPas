/**
 * @file CellOwnershipTest.cpp
 * @author D. Martin
 * @date 10.09.23
 */

#include "CellOwnershipTest.h"

using CellOwnershipTestingTypes =
    ::testing::Types<autopas::FullParticleCell<Molecule>, autopas::ReferenceParticleCell<Molecule>>;

TYPED_TEST_CASE_P(CellOwnershipTestTyped);

TEST(CellOwnershipTest, testOwnershipFullParticleCell) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::FullParticleCell<Molecule> cell({1., 1., 1.});

    cell.setPossibleParticleOwnerships(ownershipState);

    Molecule po({0.1, 0.1, 0.1}, {0., 0., 0.}, 0);
    po.setOwnershipState(autopas::OwnershipState::owned);

    Molecule ph({0.2, 0.2, 0.2}, {0., 0., 0.}, 0);
    ph.setOwnershipState(autopas::OwnershipState::halo);

    Molecule pd({0.3, 0.3, 0.3}, {0., 0., 0.}, 0);
    pd.setOwnershipState(autopas::OwnershipState::dummy);

    Molecule p1, p2;

    if (ownershipState == autopas::OwnershipState::owned) {
      p1 = po;
      p2 = ph;
    } else {
      p1 = ph;
      p2 = po;
    }

    EXPECT_NO_THROW(cell.addParticle(p1));

    // A dummy can be added always
    EXPECT_NO_THROW(cell.addParticle(pd));

    // A halo particle can not be added to a pure owned cell and vice versa
    EXPECT_ANY_THROW(cell.addParticle(p2));
  }
}

TEST(CellOwnershipTest, testOwnershipReferenceParticleCell) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::ReferenceParticleCell<Molecule> cell({1., 1., 1.});

    cell.setPossibleParticleOwnerships(ownershipState);

    Molecule po({0.1, 0.1, 0.1}, {0., 0., 0.}, 0);
    po.setOwnershipState(autopas::OwnershipState::owned);

    Molecule ph({0.2, 0.2, 0.2}, {0., 0., 0.}, 0);
    ph.setOwnershipState(autopas::OwnershipState::halo);

    Molecule pd({0.3, 0.3, 0.3}, {0., 0., 0.}, 0);
    pd.setOwnershipState(autopas::OwnershipState::dummy);

    Molecule p1, p2;

    if (ownershipState == autopas::OwnershipState::owned) {
      p1 = po;
      p2 = ph;
    } else {
      p1 = ph;
      p2 = po;
    }

    EXPECT_NO_THROW(cell.addParticleReference(&p1));

    // A dummy can be added always
    EXPECT_NO_THROW(cell.addParticleReference(&pd));

    // A halo particle can not be added to a pure owned cell and vice versa
    EXPECT_ANY_THROW(cell.addParticleReference(&p2));
  }
}

TEST(CellOwnershipTest, testOwnershipSortedCellView) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::FullParticleCell<Molecule> cell({1., 1., 1.});

    cell.setPossibleParticleOwnerships(ownershipState);

    auto scv = autopas::SortedCellView<Molecule, autopas::FullParticleCell<Molecule>>(cell, {0., 0., 0.});

    EXPECT_TRUE(scv.getPossibleParticleOwnerships() == ownershipState);

    EXPECT_ANY_THROW(scv.setPossibleParticleOwnerships(ownershipState));
  }
}

TEST(CellOwnershipTest, testOwnershipOctreeLeafNode) {
  for (auto ownershipState : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    autopas::OctreeLeafNode<Molecule> cell({0, 0, 0}, {1., 1., 1.}, nullptr, 5, 1., 1.);

    cell.setPossibleParticleOwnerships(ownershipState);

    Molecule po({0.1, 0.1, 0.1}, {0., 0., 0.}, 0);
    po.setOwnershipState(autopas::OwnershipState::owned);

    Molecule ph({0.2, 0.2, 0.2}, {0., 0., 0.}, 0);
    ph.setOwnershipState(autopas::OwnershipState::halo);

    Molecule pd({0.3, 0.3, 0.3}, {0., 0., 0.}, 0);
    pd.setOwnershipState(autopas::OwnershipState::dummy);

    Molecule p1, p2;

    if (ownershipState == autopas::OwnershipState::owned) {
      p1 = po;
      p2 = ph;
    } else {
      p1 = ph;
      p2 = po;
    }

    EXPECT_NO_THROW(cell.insert(p1));

    // A dummy can be added always
    EXPECT_NO_THROW(cell.insert(pd));

    // A halo particle can not be added to a pure owned cell and vice versa
    EXPECT_ANY_THROW(cell.insert(p2));
  }
}

TYPED_TEST_P(CellOwnershipTestTyped, testSetOwnershipState) {
  TypeParam cell({1., 1., 1.});

  // expect OwnershipState ownedOrHalo
  EXPECT_TRUE(cell.getPossibleParticleOwnerships() == (autopas::OwnershipState::owned | autopas::OwnershipState::halo));

  cell.setPossibleParticleOwnerships(autopas::OwnershipState::owned);

  // cell ownership can not be changed once set
  EXPECT_ANY_THROW(cell.setPossibleParticleOwnerships(autopas::OwnershipState::halo));
}

REGISTER_TYPED_TEST_CASE_P(CellOwnershipTestTyped, testSetOwnershipState);
INSTANTIATE_TYPED_TEST_CASE_P(TypedTest, CellOwnershipTestTyped, CellOwnershipTestingTypes);