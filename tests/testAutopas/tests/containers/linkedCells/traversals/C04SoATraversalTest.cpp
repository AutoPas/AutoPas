/**
 * @file C04SoATraversalTest.cpp
 * @author C. Menges
 * @date 06.07.2019
 */

#include "C04SoATraversalTest.h"

#include "autopas/containers/linkedCells/traversals/C04SoATraversal.h"
#include "autopas/containers/linkedCells/traversals/C08Traversal.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;

/**
 * Tests whether the interaction between two cells is correct
 */
TEST_F(C04SoATraversalTest, testTraversal) {
  std::array<size_t, 3> edgeLength = {6, 6, 6};
  autopas::LJFunctor<Molecule, FMCell> functor(1.);
  functor.setParticleProperties(24, 1);
  std::vector<FMCell> cells1;
  std::vector<FMCell> cells2;
  cells1.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);
  cells2.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  Molecule defaultParticle;
  autopasTools::generators::GridGenerator::fillWithParticles(cells1, edgeLength, {3ul, 3ul, 3ul}, defaultParticle,
                                                             {0.5, 0.5, 0.5}, {0.75, 0.75, 0.75});
  autopasTools::generators::GridGenerator::fillWithParticles(cells2, edgeLength, {3ul, 3ul, 3ul}, defaultParticle,
                                                             {0.5, 0.5, 0.5}, {0.75, 0.75, 0.75});

  NumThreadGuard numThreadGuard(1);

  autopas::C04SoATraversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, autopas::DataLayoutOption::soa, true>
      c04SoATraversal(edgeLength, &functor, 2, {1., 1., 1.});
  c04SoATraversal.setCellsToTraverse(cells1);
  c04SoATraversal.initTraversal();
  c04SoATraversal.traverseParticlePairs();
  c04SoATraversal.endTraversal();

  autopas::C08Traversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, autopas::DataLayoutOption::soa, true>
      c08Traversal(edgeLength, &functor, 2, {1., 1., 1.});
  c08Traversal.setCellsToTraverse(cells2);
  c08Traversal.initTraversal();
  c08Traversal.traverseParticlePairs();
  c08Traversal.endTraversal();

  auto iter1 = cells1.begin();
  auto iter2 = cells2.begin();
  for (; iter1 != cells1.end(); ++iter1, ++iter2) {
    if (iter1->numParticles() == 0) {
      continue;
    }
    for (auto d = 0ul; d < 3; ++d) {
      auto pos = iter1->begin()->getR();
      EXPECT_DOUBLE_EQ(iter1->begin()->getF()[d], iter2->begin()->getF()[d])
          << "Dim: " << d << " Pos: " << pos[0] << ", " << pos[1] << ", " << pos[2]
          << " ID: " << iter1->begin()->getID();
    }
  }
}
