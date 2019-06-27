/**
 * @file VerletListsCellsTraversalTest.cpp
 * @author nguyen
 * @date 26.09.18
 */

#include "VerletListsCellsTraversalTest.h"
#include "testingHelpers/NumThreadGuard.h"

VerletListsCellsTraversalTest::VerletListsCellsTraversalTest()
    : _verletListsCells(getBoxMin(), getBoxMax(), getCutoff(), autopas::TraversalOption::c18, 0.1 * getCutoff()),
      _verletListsCells_cs2(getBoxMin(), getBoxMax(), getCutoff(), autopas::TraversalOption::c18, 0.1 * getCutoff(),
                            2.0) {
  double eps = 1.0;
  double sig = 1.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
}

std::vector<unsigned long> getKernelCallsAllTraversals(autopas::VerletListsCells<Molecule> &verletListsCells,
                                                       double cutoff) {
  auto dim = verletListsCells.getCellsPerDimension();

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsC01(cutoff), flopsC18(cutoff), flopsSli(cutoff);
  autopas::FlopCounterFunctor<Molecule, FMCell> flopsC18N3(cutoff), flopsSliN3(cutoff);

  autopas::C01TraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::aos, false>
      traversalC01FLOPS(dim, &flopsC01);
  autopas::C18TraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::aos, false>
      traversalC18FLOPS(dim, &flopsC18);
  autopas::SlicedTraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::aos, false>
      traversalSliFLOPS(dim, &flopsSli);
  autopas::C18TraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::aos, true>
      traversalC18N3FLOPS(dim, &flopsC18N3);
  autopas::SlicedTraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::aos, true>
      traversalSliN3FLOPS(dim, &flopsSliN3);
  verletListsCells.rebuildNeighborLists(&traversalC01FLOPS);
  verletListsCells.iteratePairwise(&flopsC01, &traversalC01FLOPS);
  verletListsCells.rebuildNeighborLists(&traversalC18FLOPS);
  verletListsCells.iteratePairwise(&flopsC18, &traversalC18FLOPS);
  verletListsCells.rebuildNeighborLists(&traversalSliFLOPS);
  verletListsCells.iteratePairwise(&flopsSli, &traversalSliFLOPS);
  verletListsCells.rebuildNeighborLists(&traversalC18N3FLOPS);
  verletListsCells.iteratePairwise(&flopsC18N3, &traversalC18N3FLOPS);
  verletListsCells.rebuildNeighborLists(&traversalSliN3FLOPS);
  verletListsCells.iteratePairwise(&flopsSliN3, &traversalSliN3FLOPS);

  EXPECT_EQ(flopsC18.getKernelCalls(), flopsC01.getKernelCalls());
  EXPECT_EQ(flopsC18.getKernelCalls(), flopsSli.getKernelCalls());
  EXPECT_EQ(flopsC18.getKernelCalls(), (flopsC18N3.getKernelCalls() * 2));
  EXPECT_EQ(flopsC18N3.getKernelCalls(), flopsSliN3.getKernelCalls());

  return std::vector<unsigned long>({flopsC01.getKernelCalls(), flopsC18.getKernelCalls(), flopsSli.getKernelCalls(),
                                     flopsC18N3.getKernelCalls(), flopsSliN3.getKernelCalls()});
}
/**
 * Generate a VerletListCells Container and test
 * if different traversals generate the same number
 * of kernel calls.
 * @param numMolecules number of molecules in the container
 */
void VerletListsCellsTraversalTest::test(unsigned long numMolecules) {
  NumThreadGuard numThreadGuard(1);

  RandomGenerator::fillWithParticles(_verletListsCells, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0),
                                     numMolecules);
  RandomGenerator::fillWithParticles(_verletListsCells_cs2, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0),
                                     numMolecules);

  auto kCVerlet = getKernelCallsAllTraversals(_verletListsCells, getCutoff());
  auto kCVerlet_cs2 = getKernelCallsAllTraversals(_verletListsCells_cs2, getCutoff());

  EXPECT_EQ(kCVerlet, kCVerlet_cs2);
}

TEST_F(VerletListsCellsTraversalTest, test100) {
  unsigned long numMolecules = 100;
  test(numMolecules);
}

TEST_F(VerletListsCellsTraversalTest, test1000) {
  unsigned long numMolecules = 1000;
  test(numMolecules);
}

TEST_F(VerletListsCellsTraversalTest, test2000) {
  unsigned long numMolecules = 2000;
  test(numMolecules);
}
