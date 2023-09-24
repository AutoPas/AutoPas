/**
 * @file SingleCellForEachTest.cpp
 * @author lgaertner
 * @date 25.08.2021
 */

#include "SingleCellForEachTest.h"

#include <gtest/gtest.h>
#include <tests/autopasInterface/AutoPasInterfaceTest.h>

using namespace autopas;

void SingleCellForEachTest::SetUp() {
  for (int i = 0; i < 4; ++i) {
    std::array<double, 3> arr{};
    for (auto &a : arr) {
      a = static_cast<double>(i);
    }
    Molecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), 0);
    _vecOfMolecules.push_back(m);
  }

  _vecOfMolecules[0].setOwnershipState(OwnershipState::dummy);
  _vecOfMolecules[1].setOwnershipState(OwnershipState::halo);
}

void SingleCellForEachTest::TearDown() {}

/**
 * Testing logic for SingleCellForEachTest by adding particle indices to a list of 'foundParticles' and comparing to
 * 'expectedIndices'.
 * @tparam Cell FullParticleCell or ReferenceParticleCell
 * @param cell
 * @param expectedIndices of particles that should be taken into account for forEach
 * @param iteratorBehavior preferred IteratorBehavior to test cell forEach against
 * @param lowerCorner preferred lower corner of bounding box to test cell against forEachInRegion
 * @param higherCorner preferred higher corner of bounding box to test cell against forEachInRegion
 */
template <typename Cell>
void SingleCellForEachTest::testCell(Cell cell, std::vector<size_t> &expectedIndices, IteratorBehavior iteratorBehavior,
                                     std::array<double, 3> const lowerCorner,
                                     std::array<double, 3> const higherCorner) {
  std::vector<size_t> foundParticles;

  auto forEachLambda = [&](auto &p) {
    // compare the position of particle iterated over with forEach to particle with same index
    auto referenceMolecule = getMoleculeWithId(p.getID());
    for (int d = 0; d < 3; ++d) {
      EXPECT_NEAR(p.getR()[d], referenceMolecule.getR()[d], 1e-12);
    }

    // add particle to foundParticles list
    foundParticles.push_back(p.getID());
  };

  // call correct forEach function of cell (with or without region check)
  if (lowerCorner[0] == 0.f) {
    cell.forEach(forEachLambda, iteratorBehavior);
  } else {
    cell.forEach(forEachLambda, lowerCorner, higherCorner, iteratorBehavior);
  }

  ASSERT_THAT(foundParticles, ::testing::UnorderedElementsAreArray(expectedIndices));
}

/**
 * Test forEach for all particles in FullParticleCell.
 */
TEST_F(SingleCellForEachTest, testAllParticlesFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, dummy, dummy);
}

/**
 * Test forEach for owned or halo particles in FullParticleCell.
 */
TEST_F(SingleCellForEachTest, testOwnedOrHaloFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHalo, dummy, dummy);
}

/**
 * Test forEach for owned particles in FullParticleCell.
 */
TEST_F(SingleCellForEachTest, testOwnedFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);

  testCell(fpc, expectedIndices, IteratorBehavior::owned, dummy, dummy);
}

/**
 * Test forEach for all particles in a certain region of FullParticleCell.
 */
TEST_F(SingleCellForEachTest, testAllParticlesInRegionFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  const std::array<double, 3> lowerCorner{0.5, 0.5, 0.5};
  const std::array<double, 3> higherCorner{2.5, 2.5, 2.5};

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, lowerCorner, higherCorner);
}

/**
 * Test forEach for all particles of ReferenceParticleCell.
 */
TEST_F(SingleCellForEachTest, testAllParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, dummy, dummy);
}

/**
 * Test forEach for owned or halo particles of ReferenceParticleCell.
 */
TEST_F(SingleCellForEachTest, testOwnedOrHaloParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHalo, dummy, dummy);
}

/**
 * Test forEach for owned particles of ReferenceParticleCell.
 */
TEST_F(SingleCellForEachTest, testOwnedParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);

  testCell(rpc, expectedIndices, IteratorBehavior::owned, dummy, dummy);
}

/**
 * Test forEach for all particles in a certain region of ReferenceParticleCell.
 */
TEST_F(SingleCellForEachTest, testAllParticlesInRegionRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  const std::array<double, 3> lowerCorner{0.5, 0.5, 0.5};
  const std::array<double, 3> higherCorner{2.5, 2.5, 2.5};

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, lowerCorner, higherCorner);
}
