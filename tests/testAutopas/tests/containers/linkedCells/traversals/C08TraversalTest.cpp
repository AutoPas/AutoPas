/**
 * @file C08TraversalTest.cpp
 * @author F. Gratl
 * @date 24.05.18
 */

#include "C08TraversalTest.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;
using ::testing::AtLeast;

void testC08Traversal(const std::array<size_t, 3>& edgeLength, unsigned long overlap = 1ul) {
  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);
  GridGenerator::fillWithParticles<autopas::Particle>(cells, edgeLength);

  NumThreadGuard(4);

  autopas::C08Traversal<FPCell, MFunctor, false, true> c08Traversal(edgeLength, &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  const auto overlap2 = 2 * overlap;
  const long interactions = ((overlap2 + 1) * (overlap2 + 1) * (overlap2 + 1) - 1) / 2;
  EXPECT_CALL(functor, AoSFunctor(_, _, true))
      .Times((edgeLength[0] - overlap) * (edgeLength[1] - overlap) * (edgeLength[2] - overlap) * interactions);

  c08Traversal.traverseCellPairs(cells);
}

TEST_F(C08TraversalTest, testTraversal10x10x10) {
  std::array<size_t, 3> edgeLength = {10, 10, 10};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, DISABLED_testTraversal3x3x3_overlap2) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};
  testC08Traversal(edgeLength, 2);
}

TEST_F(C08TraversalTest, testTraversal2x2x2) {
  std::array<size_t, 3> edgeLength = {2, 2, 2};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, testTraversal3x3x3) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, testTraversal2x3x4) {
  std::array<size_t, 3> edgeLength = {2, 3, 4};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, testTraversal7x8x9) {
  std::array<size_t, 3> edgeLength = {7, 8, 9};
  testC08Traversal(edgeLength);
}