/**
 * @file C04SoATraversalTest.cpp
 * @author C. Menges
 * @date 06.07.2019
 */

#include "C04SoATraversalTest.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;

/**
 * Tests whether the interaction between two cells is correct
 */
TEST_F(C04SoATraversalTest, testTraversal) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};
  autopas::LJFunctor<Molecule, FMCell> functor(1.);
  functor.setParticleProperties(24, 1);
  std::vector<FMCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  Molecule default1({1.75, 2.1, 1.75}, {0., 0., 0.}, 0, 0);
  Molecule default2({1.75, 1.6, 1.75}, {0., 0., 0.}, 1, 0);

  cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength)].addParticle(default1);
  cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength)].addParticle(default2);

  NumThreadGuard numThreadGuard(1);

  autopas::C04SoATraversal<FMCell, decltype(functor), autopas::DataLayoutOption::soa, true> c04SoATraversal(
      edgeLength, &functor, 1, {1., 1., 1.});
  c04SoATraversal.setCellsToTraverse(cells);
  c04SoATraversal.initTraversal();
  c04SoATraversal.traverseParticlePairs();
  c04SoATraversal.endTraversal();

  size_t num = 0;
  for (auto &cell : cells) {
    num += cell.numParticles();
  }
  EXPECT_EQ(num, 2);
  auto iter = cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength)].begin();
  auto iter2 = cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength)].begin();
  EXPECT_DOUBLE_EQ(iter->getF()[1], -390144.0);
  EXPECT_EQ(-iter->getF()[1], iter2->getF()[1]);
}