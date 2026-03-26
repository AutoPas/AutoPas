/**
 * @file SlicedTraversalTest.cpp
 * @author F. Gratl
 * @date 01.05.18
 */

#include "SlicedTraversalTest.h"

#include "autopas/containers/linkedCells/traversals/LCSlicedTraversal.h"
#include "autopasTools/generators/GridGenerator.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;

void testSlicedTraversal(const std::array<size_t, 3> &edgeLength) {
  // Get LJ Functor with FLOP Counting enabled
  LJFunctorType</*shift*/ false, /*mixing*/ false, autopas::FunctorN3Modes::Both, /*globals*/ false, /*flops*/ true>
      ljFunctor(1.);
  std::vector<FMCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopasTools::generators::GridGenerator::fillWithParticles(cells, edgeLength, edgeLength,
                                                             autopas::utils::ParticleTypeTrait<FMCell>::value(),
                                                             {0.99, 0.99, 0.99}, {0.5, 0.5, 0.5}, {1., 1., 1.});

  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FMCell, decltype(ljFunctor)> slicedTraversal(edgeLength, &ljFunctor, 1., {1., 1., 1.},
                                                                          autopas::DataLayoutOption::aos, true);

  EXPECT_TRUE(slicedTraversal.isApplicableToDomain());
  slicedTraversal.setCellsToTraverse(cells);
  slicedTraversal.initTraversal();
  slicedTraversal.traverseParticles();

  /* center point interactions:
   *
   *           x         |         x         |        x         |                   |
   *           |  x      |         |  x      |        |  x      |            x      |
   *      x----o----x    |    x----o----x    |        o----x    |         o----x    |
   *        x  |         |         |         |        |         |         |         |
   *           x         |         x         |        x         |         x         |
   *                     |                   |                  |                   |
   *        center       |       face        |      edge        |      vertex       |
   *
   *  3x3x3 grid is a cube with 27 = 1 + 6 + 12 + 8 particles, referring to
   *  1: one center point
   *  6: six faces
   *  12: twelve edges
   *  8: eight vertices
   *
   *
   *  Center point calculates distances to 26 particles but only interacts with 6 particles
   *  "Face points" calculate distances to 17 particles but only interacts with 5 particles each
   *  "Edge points" calculate distances to 11 particles but only interacts with 4 particles each
   *  "Vertex points" calculate distances to 7 particles but only interacts with 3 particles each
   *
   *  See visualizations above for particle interaction numbers stated.
   *
   *  This leads to
   *  (1 * 26) + (6 * 17) + (12 * 11) + (8 * 7) = 316 distances calculations and
   *  (1 * 6) + (6 * 5) + (12 * 4) + (8 * 3) = 108 interactions/kernel calls.
   *  Newton3 is enabled, so there should be 316/2 = 158 distance calculations and 108/2 = 54 kernel calls (actually
   *  less, see below).
   *
   *  But: the cells in the last layer (i.e. in [2,3] for x, y, z) are considered as halo cells.
   *  This means that interactions with them have to be calculated, but not in between them.
   *  For a visualization, see https://www.geogebra.org/3d/cnuekxyk
   *
   *  - Per halo "face", there are 8 distance calculations that lead to kernel calls that are not shared with another
   *    halo "face" (see horizontal (---) and vertical ( | ) lines in visualization below)
   *  - Per halo "face", there are 8 "diagonal" distance calculations without kernel calls (see diagonal ( x ) lines in
   *    visualization below, each x = 2 lines)
   *  - => Per halo "face" there are 16 "unshared" distance calculations and 8 "unshared" kernel calls.
   *  - => 48 "unshared" distance calculations and 24 "unshared" kernel calls.
   *  - There are a further 2 distance calculations that lead to kernel calls per shared edge (see dotted lines (:, ...)
   *    in visualization below)
   *  - => 3x2 = 6 distance calculations and kernel calls
   *  - => A total of 54 distance calculations and 30 kernel calls are not calculated, as they are pure halo
   *    interactions.
   *
   *                      o...o...o
   *                      | x | x :
   *                      o---o---o
   *                      | x | x :
   *                      o---o---o
   *
   * So we expect 158-54 = 104 distance calculations, 54-30 = 24 kernel calls, leading to a hit rate of
   * 24/104=0.23076923076923078
   */
  const double expectedHitRate = 0.23076923076923078;
  ASSERT_NEAR(expectedHitRate, ljFunctor.getHitRate(), 1e-10);
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

  autopas::LCSlicedTraversal<FPCell, MPairwiseFunctor> slicedTraversal({1, 1, 1}, nullptr, 1., {1., 1., 1.},
                                                                       autopas::DataLayoutOption::aos, true);

  EXPECT_FALSE(slicedTraversal.isApplicableToDomain());
}

TEST_F(SlicedTraversalTest, testIsApplicableShrinkable) {
  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MPairwiseFunctor> slicedTraversal({5, 5, 5}, nullptr, 1., {1., 1., 1.},
                                                                       autopas::DataLayoutOption::aos, true);

  EXPECT_TRUE(slicedTraversal.isApplicableToDomain());
}

TEST_F(SlicedTraversalTest, testIsApplicableOk) {
  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MPairwiseFunctor> slicedTraversal({11, 11, 11}, nullptr, 1., {1., 1., 1.},
                                                                       autopas::DataLayoutOption::aos, true);

  EXPECT_TRUE(slicedTraversal.isApplicableToDomain());
}

TEST_F(SlicedTraversalTest, testIsApplicableOkOnlyOneDim) {
  NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, MPairwiseFunctor> slicedTraversal({1, 1, 11}, nullptr, 1., {1., 1., 1.},
                                                                       autopas::DataLayoutOption::aos, true);

  EXPECT_TRUE(slicedTraversal.isApplicableToDomain());
}
