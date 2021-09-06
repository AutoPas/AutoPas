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
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"

using ::testing::_;

void testSlicedTraversal(const std::array<size_t, 3> &edgeLength) {

  // MFunctor functor;
  autopas::FlopCounterFunctor<autopas::Particle> f(1);
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopasTools::generators::GridGenerator::fillWithParticles(cells, edgeLength, edgeLength,
                                                             autopas::utils::ParticleTypeTrait<FPCell>::value(),
                                                             {0.99, 0.99, 0.99}, {0.5, 0.5, 0.5}, {1., 1., 1.});

  //NumThreadGuard numThreadGuard(4);

  autopas::LCSlicedTraversal<FPCell, autopas::FlopCounterFunctor<autopas::Particle>, autopas::DataLayoutOption::aos, true> slicedTraversal(
      edgeLength, &f, 1., {1., 1., 1.});

  EXPECT_TRUE(slicedTraversal.isApplicable());
  slicedTraversal.setCellsToTraverse(cells);
  slicedTraversal.initTraversal();
  slicedTraversal.traverseParticlePairs();


//    EXPECT_CALL(f, FlopCounterFunctor::AoSFunctor(_, _, true))
//        .Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1) * 13);
  // all interactions: (edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1) * 13;

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
   *  Center point interacts with 6 particles (s. visualization above)
   *  "Face points" interact with 5 particles each (s. vis/ above wo/ particle a)
   *  "Edge points" interact with 4 particles each (s. vis/ above wo/ particles a and b)
   *  "Vertex points" interact with 3 particles each (s. vis/ above wo/ a, b, c)
   *  This leads to (1 * 6) + (6 * 5) + (12 * 4) + (8 * 3) = 108 interactions, since
   *  particles are in interaction length iff their position difference is a unit vector (see vis/ above).
   *  Newton3 is enabled, so there should be 108/2 = 54 kernel calls.
   */
  auto expectedKernelCalls = 54;
  ASSERT_EQ(expectedKernelCalls, f.getKernelCalls());

  /*
   * @todo: include boundary cells in LCSlicedTraversal
   * Currently, exactly those pair interactions are calculated where one of both particles is not at the "boundary" cells
   * (i.e., (0.5, 1.49, 1.49) and (0.5, 2.48, 1.49) interact, but not (0.5, 2.48, 1.49) and (0.5, 2.48, 2.48).
   * For a visualization of all interactions that currently are actually done, see https://www.geogebra.org/3d/cnuekxyk).
   * (Remark: all actual interactions were drawn there, the visualization is 100% representing the interactions.
   * The particles were drawn at the lower integer, i.e. 0.5 -> 0, 1.49 -> 1, 2.48 -> 2)
   */
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
