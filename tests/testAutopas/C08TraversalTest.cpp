/**
 * C08TraversalTest.cpp
 *
 *  Created on: 5/24/18
 *     Aauthor: F. Gratl
 */

#include "C08TraversalTest.h"

using ::testing::_;
using ::testing::AtLeast;

TEST_F(C08TraversalTest, testTraversalCube) {
  size_t edgeLength = 10;

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength * edgeLength * edgeLength);
  autopas::Particle defaultParticle;

  GridGenerator::fillWithParticles(cells, {edgeLength, edgeLength, edgeLength}, defaultParticle);
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::C08Traversal<FPCell, MFunctor, false, true> c08Traversal({edgeLength, edgeLength, edgeLength}, &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true)).Times((edgeLength - 1) * (edgeLength - 1) * (edgeLength - 1) * 13);
  c08Traversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C08TraversalTest, testTraversal2x2x2) {
  size_t edgeLength = 2;

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength * edgeLength * edgeLength);
  autopas::Particle defaultParticle;

  GridGenerator::fillWithParticles<autopas::Particle>(cells, {edgeLength, edgeLength, edgeLength}, defaultParticle);
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::C08Traversal<FPCell, MFunctor, false, true> c08Traversal({edgeLength, edgeLength, edgeLength}, &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true)).Times((edgeLength - 1) * (edgeLength - 1) * (edgeLength - 1) * 13);
  c08Traversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C08TraversalTest, testTraversal2x3x4) {
  std::array<size_t, 3> edgeLength = {2, 3, 4};

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);
  autopas::Particle defaultParticle;

  GridGenerator::fillWithParticles<autopas::Particle>(cells, {edgeLength[0], edgeLength[1], edgeLength[2]},
                                                      defaultParticle);
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