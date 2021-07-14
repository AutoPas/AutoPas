/**
 * @file SingleCellIteratorTest.cpp
 * @author tchipev
 * @date 19.01.18
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

TEST_F(SingleCellForEachTest, testAllParticlesFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy);
}

TEST_F(SingleCellForEachTest, testOwnedOrHaloFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHalo);
}

TEST_F(SingleCellForEachTest, testOwnedFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);

  testCell(fpc, expectedIndices, IteratorBehavior::owned);
}


TEST_F(SingleCellForEachTest, testAllParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy);
}

TEST_F(SingleCellForEachTest, testOwnedOrHaloParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHalo);
}

TEST_F(SingleCellForEachTest, testOwnedParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);

  testCell(rpc, expectedIndices, IteratorBehavior::owned);
}

template <typename Cell>
void SingleCellForEachTest::testCell(Cell cell, std::vector<size_t> &expectedIndices, IteratorBehavior iteratorBehavior){
  std::vector<size_t> foundParticles;

  auto forEachLambda = [&](auto &p) {
    auto referenceMolecule = getMoleculeWithId(p.getID());
    for (int d = 0; d < 3; ++d) {
      EXPECT_NEAR(p.getR()[d], referenceMolecule.getR()[d], 1e-12);
    }

    foundParticles.push_back(p.getID());
  };

  cell.forEach(forEachLambda, iteratorBehavior);

  ASSERT_THAT(foundParticles, ::testing::UnorderedElementsAreArray(expectedIndices));
}
