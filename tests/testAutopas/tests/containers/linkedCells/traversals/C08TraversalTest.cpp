/**
 * @file C08TraversalTest.cpp
 * @author F. Gratl
 * @date 24.05.18
 */

#include "C08TraversalTest.h"

using ::testing::_;
using ::testing::AtLeast;

void testC08Traversal(const std::array<size_t, 3>& edgeLength) {
  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);
  autopas::Particle defaultParticle;

  GridGenerator::fillWithParticles(cells, {edgeLength[0], edgeLength[1], edgeLength[2]}, defaultParticle);
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::C08Traversal<FPCell, MFunctor, false, true> c08Traversal({edgeLength[0], edgeLength[1], edgeLength[2]},
                                                                    &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true))
      .Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1) * 13);
  c08Traversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C08TraversalTest, testTraversal10x10x10) {
  std::array<size_t, 3> edgeLength = {10, 10, 10};
  testC08Traversal(edgeLength);
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