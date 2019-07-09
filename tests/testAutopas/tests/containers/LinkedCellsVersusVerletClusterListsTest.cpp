/**
 * @file LinkedCellsVersusVerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "LinkedCellsVersusVerletClusterListsTest.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversal.h"

template <autopas::DataLayoutOption dataLayout, bool useNewton3>
void LinkedCellsVersusVerletClusterListsTest::test(unsigned long numMolecules, double rel_err_tolerance) {
  Verlet _verletLists{getBoxMin(), getBoxMax(), getCutoff(), 0.1 * getCutoff()};
  Linked _linkedCells{getBoxMin(), getBoxMax(), getCutoff(), 0.1 * getCutoff(), 1. /*cell size factor*/};

  RandomGenerator::fillWithParticles(_linkedCells, autopas::MoleculeLJ({0., 0., 0.}, {0., 0., 0.}, 0), numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new particles
  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    _verletLists.addParticle(*it);
  }

  double eps = 1.0;
  double sig = 1.0;
  double shift = 0.0;
  autopas::MoleculeLJ::setEpsilon(eps);
  autopas::MoleculeLJ::setSigma(sig);
  autopas::LJFunctor<Molecule, FMCell> func(getCutoff(), eps, sig, shift);

  autopas::VerletClustersTraversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, dataLayout, useNewton3>
      verletTraversal(&func);
  autopas::C08Traversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, dataLayout, useNewton3> traversalLinkedLJ(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &func);
  _verletLists.rebuildNeighborLists(&verletTraversal);
  _verletLists.iteratePairwise(&verletTraversal);
  _linkedCells.iteratePairwise(&traversalLinkedLJ);

  std::vector<std::array<double, 3>> forcesVerlet(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, skip id=0 to avoid dummy particles
  for (auto it = _verletLists.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    if (m.getID() != 0) forcesVerlet.at(m.getID()) = m.getF();
  }
  forcesVerlet.at(0) = {1.0, 1.0, 1.0};

  for (auto it = _linkedCells.begin(); it.isValid(); ++it) {
    autopas::MoleculeLJ &m = *it;
    if (m.getID() != 0) forcesLinked.at(m.getID()) = m.getF();
  }
  forcesLinked.at(0) = {1.0, 1.0, 1.0};

  for (unsigned long i = 0; i < numMolecules; ++i) {
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesVerlet[i][d];
      double f2 = forcesLinked[i][d];
      double abs_err = std::abs(f1 - f2);
      double rel_err = std::abs(abs_err / f1);
      EXPECT_LT(rel_err, rel_err_tolerance);
    }
  }

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsVerlet(getCutoff()), flopsLinked(getCutoff());
  autopas::C08Traversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, dataLayout, useNewton3> traversalFLOPS(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &flopsLinked);
  autopas::VerletClustersTraversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, dataLayout, useNewton3>
      traversalFLOPSVerlet(&flopsVerlet);
  _verletLists.iteratePairwise(&traversalFLOPSVerlet);
  _linkedCells.iteratePairwise(&traversalFLOPS);

  if (not(dataLayout == autopas::DataLayoutOption::soa && not useNewton3)) {
    ASSERT_EQ(flopsLinked.getKernelCalls(), flopsVerlet.getKernelCalls());
  }
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-14;

  test<autopas::DataLayoutOption::aos, false>(numMolecules, rel_err_tolerance);
  test<autopas::DataLayoutOption::soa, false>(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 2e-12;
  test<autopas::DataLayoutOption::aos, false>(numMolecules, rel_err_tolerance);
  test<autopas::DataLayoutOption::soa, false>(numMolecules, rel_err_tolerance);
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, test2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;
  test<autopas::DataLayoutOption::aos, false>(numMolecules, rel_err_tolerance);
  test<autopas::DataLayoutOption::soa, false>(numMolecules, rel_err_tolerance);
}
