/**
 * @file C04SoATraversalTest.cpp
 * @author C. Menges
 * @date 06.07.2019
 */

#include "C04SoATraversalTest.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;

void testC04SoATraversal(const std::array<size_t, 3> &edgeLength, unsigned long overlap = 1ul) {
  autopas::LJFunctor<autopas::Particle, FPCell> functor(1., 1., 1., 1.);
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopas::Particle defaultParticle;
  GridGenerator::fillWithParticles(cells, edgeLength, defaultParticle, {1., 1., 1.}, {0.25, 0.25, 0.25});
  GridGenerator::fillWithParticles(cells, edgeLength, defaultParticle, {1., 1., 1.}, {0.75, 0.75, 0.75});

  NumThreadGuard numThreadGuard(1);

  autopas::C04SoATraversal<FPCell, autopas::LJFunctor<autopas::Particle, FPCell>, autopas::DataLayoutOption::soa, true>
      c04SoATraversal(edgeLength, &functor, overlap);
  c04SoATraversal.initTraversal(cells);
  c04SoATraversal.traverseCellPairs(cells);
  c04SoATraversal.endTraversal(cells);
  size_t innerCell1 = autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength);
  size_t innerCell2 = autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength);
  // Forces inside the domain should be equivalent
  const auto refForce = cells[innerCell1].begin()->getF();
  for (auto cellIndex : {innerCell1, innerCell2}) {
    auto iter = cells[cellIndex].begin();
    for (; iter.isValid(); ++iter) {
      for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(refForce[d], iter->getF()[d]) << "Dim: " << d << " CellIndex: " << cellIndex;
      }
    }
  }
}

TEST_F(C04SoATraversalTest, testTraversal) {
  std::array<size_t, 3> edgeLength = {3, 4, 3};
  testC04SoATraversal(edgeLength);
}