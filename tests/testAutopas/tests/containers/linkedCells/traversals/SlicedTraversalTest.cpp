/**
 * @file SlicedTraversalTest.cpp
 * @author F. Gratl
 * @date 01.05.18
 */

#include "SlicedTraversalTest.h"

#include "autopas/containers/linkedCells/traversals/LCSlicedTraversal.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;

void testSlicedTraversal(const std::array<size_t, 3> &edgeLength) {
  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopasTools::generators::GridGenerator::fillWithParticles(cells, edgeLength, edgeLength);

  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal(
      edgeLength, &functor, 1., {1., 1., 1.});

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_TRUE(slicedTraversal.isApplicable());
  EXPECT_CALL(functor, AoSFunctor(_, _, true))
      .Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1) * 13);
  slicedTraversal.setCellsToTraverse(cells);
  slicedTraversal.initTraversal();
  slicedTraversal.traverseParticlePairs();
}

/**
 * This test is the same as testTraversalCube except that the domain is too small for 4 threads.
 * It expects the sliced traversal to start less threads but still work.
 */
TEST_F(SlicedTraversalTest, testTraversalCubeShrink) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};
  testSlicedTraversal(edgeLength);
}

TEST_F(SlicedTraversalTest, testIsApplicableTooSmall) {
  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal({1, 1, 1}, nullptr,
                                                                                                     1., {1., 1., 1.});

  EXPECT_FALSE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableShrinkable) {
  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal({5, 5, 5}, nullptr,
                                                                                                     1., {1., 1., 1.});

  EXPECT_TRUE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableOk) {
  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal(
      {11, 11, 11}, nullptr, 1., {1., 1., 1.});

  EXPECT_TRUE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableOkOnlyOneDim) {
  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal(
      {1, 1, 11}, nullptr, 1., {1., 1., 1.});

  EXPECT_TRUE(slicedTraversal.isApplicable());
}
