/**
 * @file C01TraversalTest.cpp
 * @author S. Seckler
 * @date 10.01.2019
 */

#include "C01TraversalTest.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;
using ::testing::AtLeast;

void testC01Traversal(const std::array<size_t, 3>& edgeLength, unsigned long interactions, double cutoff = 1.0) {
  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  GridGenerator::fillWithParticles<autopas::Particle>(cells, edgeLength);

  NumThreadGuard(4);

  autopas::C01Traversal<FPCell, MFunctor, false, false> C01Traversal(edgeLength, &functor, cutoff);

  EXPECT_CALL(functor, AoSFunctor(_, _, false)).Times(interactions);
  C01Traversal.traverseCellPairs(cells);
}

unsigned long getNumInteractions(const std::array<size_t, 3>& edgeLength, unsigned long overlap = 1ul) {
  // every particle interacts with 26 others. First and last layer of each dim is covered by previous interactions
  const auto overlap2 = 2 * overlap;
  const auto interactionsPerCell = (overlap2 + 1) * (overlap2 + 1) * (overlap2 + 1) - 1;
  return (edgeLength[0] - overlap2) * (edgeLength[1] - overlap2) * (edgeLength[2] - overlap2) * interactionsPerCell;
}

TEST_F(C01TraversalTest, testTraversal10x10x10) {
  std::array<size_t, 3> edgeLength = {10, 10, 10};
  testC01Traversal(edgeLength, getNumInteractions(edgeLength));
}

TEST_F(C01TraversalTest, testTraversal10x10x10_overlap2) {
  std::array<size_t, 3> edgeLength = {10, 10, 10};
  testC01Traversal(edgeLength, getNumInteractions(edgeLength, 2ul), 2.0);
}

TEST_F(C01TraversalTest, testTraversal10x10x10_overlap5) {
  std::array<size_t, 3> edgeLength = {11, 11, 11};
  // number of interacting cells per cell, without cell itself: ((2*5)+1)^3 - 1 = 1330
  // number of interactions removed due to cutoff: 244
  // 1330 - 244 = 1086
  // number of calculated cells (11-2*5)^3 = 1
  testC01Traversal(edgeLength, 1 * 1086, 5.0);
}

TEST_F(C01TraversalTest, testTraversal2x2x2) {
  std::array<size_t, 3> edgeLength = {2, 2, 2};
  testC01Traversal(edgeLength, getNumInteractions(edgeLength));
}

TEST_F(C01TraversalTest, testTraversal3x3x3) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};
  testC01Traversal(edgeLength, getNumInteractions(edgeLength));
}

TEST_F(C01TraversalTest, testTraversal2x3x4) {
  std::array<size_t, 3> edgeLength = {2, 3, 4};
  testC01Traversal(edgeLength, getNumInteractions(edgeLength));
}

TEST_F(C01TraversalTest, testTraversal7x8x9) {
  std::array<size_t, 3> edgeLength = {7, 8, 9};
  testC01Traversal(edgeLength, getNumInteractions(edgeLength));
}

TEST_F(C01TraversalTest, testTraversal7x8x9_overlap2) {
  std::array<size_t, 3> edgeLength = {7, 8, 9};
  testC01Traversal(edgeLength, getNumInteractions(edgeLength, 2ul), 2.0);
}

TEST_F(C01TraversalTest, testTraversal7x8x9_overlap3) {
  std::array<size_t, 3> edgeLength = {7, 8, 9};
  // number of interacting cells per cell, without cell itself: ((2*3)+1)^3 - 1 = 342
  // number of interactions removed due to cutoff: 8
  // 342 - 8 = 334
  // number of calculated cells (7-2*3) * (8-2*3) * (9-2*3) = 1 * 2 * 3 = 6
  testC01Traversal(edgeLength, 334 * 6, 3.0);
}