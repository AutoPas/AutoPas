/**
 * @file SlicedTraversalTest.cpp
 * @author F. Gratl
 * @date 01.05.18
 */

#include "SlicedTraversalTest.h"
#include "testingHelpers/GridGenerator.h"

using ::testing::_;
using ::testing::AtLeast;

void testSlicedTraversal(const std::array<size_t, 3>& edgeLength) {
  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopas::Particle defaultParticle;
  GridGenerator::fillWithParticles(cells, edgeLength, defaultParticle);
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal(
      {edgeLength[0], edgeLength[1], edgeLength[2]}, &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true))
      .Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1) * 13);
  slicedTraversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(SlicedTraversalTest, testTraversalCube) {
  std::array<size_t, 3> edgeLength = {10, 10, 10};
  testSlicedTraversal(edgeLength);
}

/**
 * This test is the same as testTraversalCube except that the domain is too small for 4 threads.
 * It expects the sliced traversal to start less threads but still work.
 */
TEST_F(SlicedTraversalTest, testTraversalCubeShrink) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};
  testSlicedTraversal(edgeLength);
}

TEST_F(SlicedTraversalTest, testTraversalCuboid) {
  std::array<size_t, 3> edgeLength = {5, 7, 10};
  testSlicedTraversal(edgeLength);
}

TEST_F(SlicedTraversalTest, testIsApplicableTooSmall) {
  std::vector<FPCell> cells;

#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal({1, 1, 1}, nullptr);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_FALSE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableShrinkable) {
  std::vector<FPCell> cells;

#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal({5, 5, 5}, nullptr);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_TRUE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableOk) {
  std::vector<FPCell> cells;

#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal({11, 11, 11},
                                                                                                   nullptr);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_TRUE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableOkOnlyOneDim) {
  std::vector<FPCell> cells;

#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> slicedTraversal({1, 1, 11}, nullptr);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_TRUE(slicedTraversal.isApplicable());
}
