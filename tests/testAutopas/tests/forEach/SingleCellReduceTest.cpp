/**
 * @file SingleCellReduceTest.cpp
 * @author lgaertner
 * @date 25.08.2021
 */

#include "SingleCellReduceTest.h"

using namespace autopas;

void SingleCellReduceTest::SetUp() {
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

void SingleCellReduceTest::TearDown() {}

/**
 * Testing logic for SingleCellReduceTest which sums up the indices of all found particles and compares it with the sum
 * of the given expectedIndices.
 * @tparam Cell FullParticleCell or ReferenceParticleCell
 * @param cell
 * @param expectedIndices of particles that should be taken into account for reduction
 * @param iteratorBehavior preferred IteratorBehavior to test cell reduction against
 * @param lowerCorner preferred lower corner of bounding box to test cell against reduceInRegion
 * @param higherCorner preferred higher corner of bounding box to test cell against reduceInRegion
 */
template <typename Cell>
void SingleCellReduceTest::testCell(Cell cell, std::vector<size_t> &expectedIndices, IteratorBehavior iteratorBehavior,
                                    std::array<double, 3> const lowerCorner, std::array<double, 3> const higherCorner) {
  // expectedReductionValue is sum of all expectedIndices
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0ul);

  std::vector<size_t> foundParticles;
  double reductionValue = 0.0f;

  auto reduceLambda = [&](auto &p, auto &result) {
    // actual reduction
    result += p.getID();

    // similar foundParticles test as in SingleCellForEach to improve debugging
    foundParticles.push_back(p.getID());
  };

  // call correct reduction function of cell (with or without region check)
  if (lowerCorner[0] == 0.f) {
    cell.reduce(reduceLambda, reductionValue, iteratorBehavior);
  } else {
    cell.reduce(reduceLambda, reductionValue, lowerCorner, higherCorner, iteratorBehavior);
  }

  ASSERT_THAT(foundParticles, ::testing::UnorderedElementsAreArray(expectedIndices));
  EXPECT_EQ(reductionValue, expectedReductionValue);
}

/**
 * Test reduction of all particles in a FullParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testAllParticlesFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, dummy, dummy);
}

/**
 * Test reduction of owned or halo particles in a FullParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testOwnedOrHaloFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHalo, dummy, dummy);
}

/**
 * Test reduction of owned particles in a FullParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testOwnedFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);

  testCell(fpc, expectedIndices, IteratorBehavior::owned, dummy, dummy);
}

/**
 * Test reduction of all particles in a certain region of a FullParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testAllParticlesInRegionFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  const std::array<double, 3> lowerCorner{0.5, 0.5, 0.5};
  const std::array<double, 3> higherCorner{2.5, 2.5, 2.5};

  testCell(fpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, lowerCorner, higherCorner);
}

/**
 * Test reduction of all particles in a ReferenceParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testAllParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, dummy, dummy);
}

/**
 * Test reduction of owned or halo particles in a ReferenceParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testOwnedOrHaloParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHalo, dummy, dummy);
}

/**
 * Test reduction of owned particles in a ReferenceParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testOwnedParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);

  testCell(rpc, expectedIndices, IteratorBehavior::owned, dummy, dummy);
}

/**
 * Test reduction of all particles in a certain region of a ReferenceParticleCell by using testCell.
 */
TEST_F(SingleCellReduceTest, testAllParticlesInRegionRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);

  const std::array<double, 3> lowerCorner{0.5, 0.5, 0.5};
  const std::array<double, 3> higherCorner{2.5, 2.5, 2.5};

  testCell(rpc, expectedIndices, IteratorBehavior::ownedOrHaloOrDummy, lowerCorner, higherCorner);
}
