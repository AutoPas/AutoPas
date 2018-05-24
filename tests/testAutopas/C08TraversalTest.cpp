/**
 * C08TraversalTest.cpp
 *
 *  Created on: 5/24/18
 *     Aauthor: F. Gratl
 */


#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif
#include "C08TraversalTest.h"

using ::testing::_;
using ::testing::AtLeast;

void C08TraversalTest::fillWithParticles(
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

TEST_F(C08TraversalTest, testTraversalCube) {

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
  autopas::C08Traversal<FPCell, MCellFunctor> c08Traversal (cells, {edgeLength,edgeLength,edgeLength}, &cellFunctor);

  // every particle interacts with 13 others. Last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _)).Times((edgeLength - 1) * (edgeLength - 1) * (edgeLength - 1)  * 13);
  c08Traversal.traverseCellPairs();
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C08TraversalTest, testTraversal2x2x2) {

  size_t edgeLength = 2;

  MFunctor functor;
  MCellFunctor cellFunctor(&functor);
  std::vector<FPCell> cells;
  cells.resize(edgeLength*edgeLength*edgeLength);

  fillWithParticles(cells, {edgeLength,edgeLength,edgeLength});
#ifdef _OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::C08Traversal<FPCell, MCellFunctor> c08Traversal (cells, {edgeLength,edgeLength,edgeLength}, &cellFunctor);

  // every particle interacts with 13 others. Last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _)).Times((edgeLength - 1) * (edgeLength - 1) * (edgeLength - 1)  * 13);
  c08Traversal.traverseCellPairs();
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C08TraversalTest, testTraversal2x3x4) {

  std::array<size_t, 3> edgeLength = {2,3,4};

  MFunctor functor;
  MCellFunctor cellFunctor(&functor);
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0]*edgeLength[1]*edgeLength[2]);

  fillWithParticles(cells, {edgeLength[0],edgeLength[1],edgeLength[2]});
#ifdef _OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::C08Traversal<FPCell, MCellFunctor> c08Traversal(cells, {edgeLength[0],edgeLength[1],edgeLength[2]}, &cellFunctor);

  // every particle interacts with 13 others. Last layer of each dim is covered by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _)).Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1)  * 13);
  c08Traversal.traverseCellPairs();
#ifdef _OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}