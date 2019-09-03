/**
 * @file C04SoATraversalTest.cpp
 * @author C. Menges
 * @date 06.07.2019
 */

#include "C04SoATraversalTest.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;

/**
 * Tests whether the interaction between two cells is correct
 */
TEST_F(C04SoATraversalTest, testTraversal) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};

  autopas::LJFunctor<autopas::Particle, FPCell> functor(1., 1., 1., 1.);
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopas::Particle defaultParticle;
  defaultParticle.setR({1.75, 2.1, 1.75});
  cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength)].addParticle(defaultParticle);
  defaultParticle.setR({1.75, 1.6, 1.75});
  cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength)].addParticle(defaultParticle);

  NumThreadGuard numThreadGuard(1);

  autopas::C04SoATraversal<FPCell, autopas::LJFunctor<autopas::Particle, FPCell>, autopas::DataLayoutOption::soa, true>
      c04SoATraversal(edgeLength, &functor, 1, {1., 1., 1.});
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