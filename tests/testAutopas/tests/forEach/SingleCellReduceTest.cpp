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

TEST_F(SingleCellReduceTest, testAllParticlesFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  testCell(fpc, expectedIndices, expectedReductionValue, IteratorBehavior::ownedOrHaloOrDummy, dummy, dummy);
}

TEST_F(SingleCellReduceTest, testOwnedOrHaloFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  testCell(fpc, expectedIndices, expectedReductionValue, IteratorBehavior::ownedOrHalo, dummy, dummy);
}

TEST_F(SingleCellReduceTest, testOwnedFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  testCell(fpc, expectedIndices, expectedReductionValue, IteratorBehavior::owned, dummy, dummy);
}

TEST_F(SingleCellReduceTest, testAllParticlesInRegionFpc) {
  FMCell fpc;
  fillWithParticles(&fpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  const std::array<double, 3> lowerCorner{0.5, 0.5, 0.5};
  const std::array<double, 3> higherCorner{2.5, 2.5, 2.5};

  testCell(fpc, expectedIndices, expectedReductionValue, IteratorBehavior::ownedOrHaloOrDummy, lowerCorner, higherCorner);
}

TEST_F(SingleCellReduceTest, testAllParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size());
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  testCell(rpc, expectedIndices, expectedReductionValue, IteratorBehavior::ownedOrHaloOrDummy, dummy, dummy);
}

TEST_F(SingleCellReduceTest, testOwnedOrHaloParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 1);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  testCell(rpc, expectedIndices, expectedReductionValue, IteratorBehavior::ownedOrHalo, dummy, dummy);
}

TEST_F(SingleCellReduceTest, testOwnedParticlesRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 2);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  testCell(rpc, expectedIndices, expectedReductionValue, IteratorBehavior::owned, dummy, dummy);
}

TEST_F(SingleCellReduceTest, testAllParticlesInRegionRpc) {
  ReferenceParticleCell<Molecule> rpc;
  fillWithParticleReferences(&rpc);

  std::vector<size_t> expectedIndices(_vecOfMolecules.size() - 2);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 1);
  size_t expectedReductionValue = std::accumulate(expectedIndices.begin(), expectedIndices.end(), 0l);

  const std::array<double, 3> lowerCorner{0.5, 0.5, 0.5};
  const std::array<double, 3> higherCorner{2.5, 2.5, 2.5};

  testCell(rpc, expectedIndices, expectedReductionValue, IteratorBehavior::ownedOrHaloOrDummy, lowerCorner, higherCorner);
}
template <typename Cell>
void SingleCellReduceTest::testCell(Cell cell, std::vector<size_t> &expectedIndices, size_t expectedReductionValue, IteratorBehavior iteratorBehavior,
                                     std::array<double, 3> const lowerCorner,
                                     std::array<double, 3> const higherCorner) {
  std::vector<size_t> foundParticles;
  double reductionValue = 0.0f;

  auto reduceLambda = [&](auto &p, auto &result) {
    result += p.getID();

    foundParticles.push_back(p.getID());
  };

  if (lowerCorner[0] == 0.f) {
    cell.reduce(reduceLambda, reductionValue, iteratorBehavior);
  } else {
    cell.reduce(reduceLambda, reductionValue, lowerCorner, higherCorner, iteratorBehavior);
  }

  ASSERT_THAT(foundParticles, ::testing::UnorderedElementsAreArray(expectedIndices));
  EXPECT_EQ(reductionValue, expectedReductionValue);
}
