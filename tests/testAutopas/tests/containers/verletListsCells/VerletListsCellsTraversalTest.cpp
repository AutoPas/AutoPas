/**
 * @file VerletListsCellsTraversalTest.cpp
 * @author nguyen
 * @date 26.09.18
 */

#include "VerletListsCellsTraversalTest.h"

VerletListsCellsTraversalTest::VerletListsCellsTraversalTest()
    : _verletListsCells(getBoxMin(), getBoxMax(), getCutoff(), autopas::TraversalOption::c18, 0.1 * getCutoff(), 2) {
  double eps = 1.0;
  double sig = 1.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
}
/**
 * Generate a VerletListCells Container and test
 * if different traversals generate the same number
 * of kernel calls.
 * @param numMolecules number of molecules in the container
 */
void VerletListsCellsTraversalTest::test(unsigned long numMolecules) {
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  RandomGenerator::fillWithParticles(_verletListsCells, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0),
                                     numMolecules);

  auto dim = _verletListsCells.getCellsPerDimension();

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsC01(getCutoff()), flopsC18(getCutoff()), flopsSli(getCutoff());
  autopas::FlopCounterFunctor<Molecule, FMCell> flopsC18N3(getCutoff()), flopsSliN3(getCutoff());
  autopas::C01TraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos,
                              false>
      traversalC01FLOPS(dim, &flopsC01);
  autopas::C18TraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos,
                              false>
      traversalC18FLOPS(dim, &flopsC18);
  autopas::SlicedTraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos,
                                 false>
      traversalSliFLOPS(dim, &flopsSli);
  autopas::C18TraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos,
                              true>
      traversalC18N3FLOPS(dim, &flopsC18N3);
  autopas::SlicedTraversalVerlet<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, autopas::DataLayoutOption::aos,
                                 true>
      traversalSliN3FLOPS(dim, &flopsSliN3);
  _verletListsCells.iteratePairwiseAoS(&flopsC01, &traversalC01FLOPS, false);
  _verletListsCells.iteratePairwiseAoS(&flopsC18, &traversalC18FLOPS, false);
  _verletListsCells.iteratePairwiseAoS(&flopsSli, &traversalSliFLOPS, false);
  _verletListsCells.iteratePairwiseAoS(&flopsC18N3, &traversalC18N3FLOPS, true);
  _verletListsCells.iteratePairwiseAoS(&flopsSliN3, &traversalSliN3FLOPS, true);

  ASSERT_EQ(flopsC18.getKernelCalls(), flopsC01.getKernelCalls());
  ASSERT_EQ(flopsC18.getKernelCalls(), flopsSli.getKernelCalls());
  ASSERT_EQ(flopsC18.getKernelCalls(), (flopsC18N3.getKernelCalls() * 2));
  ASSERT_EQ(flopsC18N3.getKernelCalls(), flopsSliN3.getKernelCalls());

#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
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
