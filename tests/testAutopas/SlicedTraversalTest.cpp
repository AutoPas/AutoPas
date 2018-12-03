/**
 * @file SlicedTraversalTest.cpp
 * @author F. Gratl
 * @date 01.05.18
 */

#include "SlicedTraversalTest.h"

using ::testing::_;
using ::testing::AtLeast;

void SlicedTraversalTest::fillWithParticles(std::vector<FPCell> &cells, std::array<size_t, 3> particlesPerDim) {
  size_t id = 0;
  size_t cellId = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        auto p = autopas::Particle({x + .5, y + .5, z + .5}, {0, 0, 0}, id++);
        cells[cellId++].addParticle(p);
      }
    }
  }
}

TEST_F(SlicedTraversalTest, testTraversalCube) {
  size_t edgeLength = 10;

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength * edgeLength * edgeLength);

  fillWithParticles(cells, {edgeLength, edgeLength, edgeLength});
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, false, true> slicedTraversal({edgeLength, edgeLength, edgeLength},
                                                                          &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true)).Times((edgeLength - 1) * (edgeLength - 1) * (edgeLength - 1) * 13);
  slicedTraversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

/**
 * This test is the same as testTraversalCube except that the domain is too small for 4 threads.
 * It expects the sliced traversal to start less threads but still work.
 */
TEST_F(SlicedTraversalTest, testTraversalCubeShrink) {
  size_t edgeLength = 3;

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength * edgeLength * edgeLength);

  fillWithParticles(cells, {edgeLength, edgeLength, edgeLength});
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, false, true> slicedTraversal({edgeLength, edgeLength, edgeLength},
                                                                          &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true)).Times((edgeLength - 1) * (edgeLength - 1) * (edgeLength - 1) * 13);
  slicedTraversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(SlicedTraversalTest, testTraversalCuboid) {
  std::array<size_t, 3> edgeLength = {5, 7, 10};

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  fillWithParticles(cells, {edgeLength[0], edgeLength[1], edgeLength[2]});
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, false, true> slicedTraversal({edgeLength[0], edgeLength[1], edgeLength[2]},
                                                                          &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true))
      .Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1) * 13);
  slicedTraversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(SlicedTraversalTest, testIsApplicableTooSmall) {
  std::vector<FPCell> cells;

#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MFunctor, false, true> slicedTraversal({1, 1, 1}, nullptr);
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
  autopas::SlicedTraversal<FPCell, MFunctor, false, true> slicedTraversal({5, 5, 5}, nullptr);
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
  autopas::SlicedTraversal<FPCell, MFunctor, false, true> slicedTraversal({11, 11, 11}, nullptr);
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
  autopas::SlicedTraversal<FPCell, MFunctor, false, true> slicedTraversal({1, 1, 11}, nullptr);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_TRUE(slicedTraversal.isApplicable());
}