/*
 * SlicedTraversalTest.cpp
 *
 *  Created on: 22 Jan 2018
 *      Author: gratl
 */
#ifdef _OPENMP
#include <omp.h>
#endif
#include "SlicedTraversalTest.h"

using ::testing::_;
using ::testing::AtLeast;

void SlicedTraversalTest::fillWithParticles(
    std::vector<FPCell> &cells,
    std::array<size_t, 3> particlesPerDim) {
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

  // TODO: more propper way to do this would be through a mock of the cellFunctor
  size_t edgeLength = 10;

  MFunctor functor;
  MCellFunctor cellFunctor(&functor);
  std::vector<FPCell> cells;
  cells.resize(edgeLength*edgeLength*edgeLength);

  fillWithParticles(cells, {edgeLength,edgeLength,edgeLength});
#ifdef _OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {edgeLength,edgeLength,edgeLength}, &cellFunctor);

  // every particle interacts with 13 others. Last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _)).Times((edgeLength - 1) * (edgeLength - 1) * (edgeLength - 1)  * 13);
  slicedTraversal.traverseCellPairs();
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(SlicedTraversalTest, testTraversalCuboid) {

  // TODO: more propper way to do this would be through a mock of the cellFunctor
  std::array<size_t, 3> edgeLength = {5,7,10};

  MFunctor functor;
  MCellFunctor cellFunctor(&functor);
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0]*edgeLength[1]*edgeLength[2]);

  fillWithParticles(cells, {edgeLength[0],edgeLength[1],edgeLength[2]});
#ifdef _OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {edgeLength[0],edgeLength[1],edgeLength[2]}, &cellFunctor);

  // every particle interacts with 13 others. Last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _)).Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1)  * 13);
  slicedTraversal.traverseCellPairs();
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(SlicedTraversalTest, testIsApplicableTooSmall) {
  std::vector<FPCell> cells;

#ifdef _OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {1,1,1}, nullptr);
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_FALSE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableOk) {
  std::vector<FPCell> cells;

#ifdef _OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {11,11,11}, nullptr);
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_TRUE(slicedTraversal.isApplicable());
}

TEST_F(SlicedTraversalTest, testIsApplicableOkOnlyOneDim) {

  std::vector<FPCell> cells;

#ifdef _OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::SlicedTraversal<FPCell, MCellFunctor> slicedTraversal(cells, {1,1,11}, nullptr);
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif

  EXPECT_TRUE(slicedTraversal.isApplicable());
}