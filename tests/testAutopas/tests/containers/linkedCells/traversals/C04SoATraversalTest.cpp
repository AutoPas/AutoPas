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
double universalValue=1; //epsilon=sigma=mass=1.0
ParticleClassLibrary PCL = ParticleClassLibrary(universalValue,universalValue,universalValue);
  autopas::LJFunctor<Molecule, FMCell> functor(1., PCL, 1.);
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  autopas::MoleculeLJ<> default1 = Molecule({1.75, 2.1, 1.75},{0.,0.,0.},0);
    autopas::MoleculeLJ<> default2 = Molecule({1.75, 1.6, 1.75},{0.,0.,0.},1);

    cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength)].addParticle(default1);
  cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength)].addParticle(default2);

  NumThreadGuard numThreadGuard(1);

  autopas::C04SoATraversal<FPCell, autopas::LJFunctor<Molecule, FMCell>, autopas::DataLayoutOption::soa, true>
      c04SoATraversal(edgeLength, &functor, 1);
  c04SoATraversal.setCellsToTraverse(cells);
  c04SoATraversal.initTraversal();
  c04SoATraversal.traverseParticlePairs();
  c04SoATraversal.endTraversal();

  size_t num = 0;
  for (auto cell : cells) {
    num += cell.numParticles();
  }
  EXPECT_EQ(num, 2);
  auto iter = cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, edgeLength)].begin();
  auto iter2 = cells[autopas::utils::ThreeDimensionalMapping::threeToOneD(1ul, 2ul, 1ul, edgeLength)].begin();
  EXPECT_DOUBLE_EQ(iter->getF()[1], -390144.0);
  EXPECT_EQ(-iter->getF()[1], iter2->getF()[1]);
}