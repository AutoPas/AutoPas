/**
 * @file C01TraversalTest.cpp
 * @author S. Seckler
 * @date 10.01.2019
 */

#include "C01TraversalTest.h"

using ::testing::_;
using ::testing::AtLeast;

TEST_F(C01TraversalTest, testTraversalCube) {
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
  autopas::C01Traversal<FPCell, MFunctor, false, false> C01Traversal({edgeLength, edgeLength, edgeLength}, &functor);

  // every particle interacts with 26 others. First and last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, false)).Times((edgeLength - 2) * (edgeLength - 2) * (edgeLength - 2) * 26);
  C01Traversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C01TraversalTest, testTraversal2x2x2) {
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
  autopas::C01Traversal<FPCell, MFunctor, false, false> C01Traversal({edgeLength, edgeLength, edgeLength}, &functor);

  // every particle interacts with 26 others. First and last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, false)).Times((edgeLength - 2) * (edgeLength - 2) * (edgeLength - 2) * 26);
  C01Traversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C01TraversalTest, testTraversal3x3x4) {
  std::array<size_t, 3> edgeLength = {3, 3, 4};

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
  autopas::C01Traversal<FPCell, MFunctor, false, false> C01Traversal({edgeLength[0], edgeLength[1], edgeLength[2]},
                                                                    &functor);

  // every particle interacts with 26 others. First and last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, false))
      .Times((edgeLength[0] - 2) * (edgeLength[1] - 2) * (edgeLength[2] - 2) * 26);
  C01Traversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}