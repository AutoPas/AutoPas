/**
 * @file LinkedCellsVersusVerletClusterListsTest.cpp
 * @author nguyen
 * @date 21.10.18
 */

#include "LinkedCellsVersusVerletClusterListsTest.h"
#include <autopas/selectors/TraversalSelector.h>
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversal.h"

LinkedCellsVersusVerletClusterListsTest::LinkedCellsVersusVerletClusterListsTest() {}

template <autopas::DataLayoutOption dataLayout, bool useNewton3>
void LinkedCellsVersusVerletClusterListsTest::test(unsigned long numMolecules, double rel_err_tolerance,
                                                   autopas::TraversalOption traversalOption,
                                                   std::array<double, 3> boxMax) {
  Verlet _verletLists{getBoxMin(), boxMax, getCutoff(), 0.1 * getCutoff(), 2};
  Linked _linkedCells{getBoxMin(), boxMax, getCutoff()};

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

  auto verletTraversal = autopas::TraversalSelector<FMCell>::generateTraversal(
      traversalOption, func, _verletLists.getTraversalSelectorInfo(), dataLayout,
      useNewton3 ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled);

  autopas::C08Traversal<FMCell, autopas::LJFunctor<Molecule, FMCell>, dataLayout, useNewton3> traversalLinkedLJ(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &func);

  _verletLists.iteratePairwise(&func, &*verletTraversal);
  _linkedCells.iteratePairwise(&func, &traversalLinkedLJ);

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
      EXPECT_NEAR(f1, f2, std::fabs(f1 * rel_err_tolerance));
    }
  }

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsVerlet(getCutoff()), flopsLinked(getCutoff());

  autopas::C08Traversal<FMCell, autopas::FlopCounterFunctor<Molecule, FMCell>, dataLayout, useNewton3> traversalFLOPS(
      _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &flopsLinked);

  auto traversalFLOPSVerlet = autopas::TraversalSelector<FMCell>::generateTraversal(
      traversalOption, flopsVerlet, _verletLists.getTraversalSelectorInfo(), dataLayout,
      useNewton3 ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled);

  _verletLists.iteratePairwise(&flopsVerlet, &*traversalFLOPSVerlet);
  _linkedCells.iteratePairwise(&flopsLinked, &traversalFLOPS);

  if (not(dataLayout == autopas::DataLayoutOption::soa && not useNewton3)) {
    ASSERT_EQ(flopsLinked.getKernelCalls(), flopsVerlet.getKernelCalls());
  }
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, verletClustersTest100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-14;

  for (auto boxMax : {getBoxMaxBig(), getBoxMaxSmall()}) {
    test<autopas::DataLayoutOption::aos, false>(numMolecules, rel_err_tolerance,
                                                autopas::TraversalOption::verletClusters, boxMax);
    test<autopas::DataLayoutOption::soa, false>(numMolecules, rel_err_tolerance,
                                                autopas::TraversalOption::verletClusters, boxMax);
  }
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, verletClustersTest1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 2e-12;

  for (auto boxMax : {getBoxMaxBig(), getBoxMaxSmall()}) {
    test<autopas::DataLayoutOption::aos, false>(numMolecules, rel_err_tolerance,
                                                autopas::TraversalOption::verletClusters, boxMax);
    test<autopas::DataLayoutOption::soa, false>(numMolecules, rel_err_tolerance,
                                                autopas::TraversalOption::verletClusters, boxMax);
  }
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, verletClustersTest2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;

  for (auto boxMax : {getBoxMaxBig(), getBoxMaxSmall()}) {
    test<autopas::DataLayoutOption::soa, false>(numMolecules, rel_err_tolerance,
                                                autopas::TraversalOption::verletClusters, boxMax);
    test<autopas::DataLayoutOption::aos, false>(numMolecules, rel_err_tolerance,
                                                autopas::TraversalOption::verletClusters, boxMax);
  }
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, verletClustersColoringTest100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1.5e-14;

  for (auto boxMax : {getBoxMaxBig(), getBoxMaxSmall()}) {
    test<autopas::DataLayoutOption::aos, true>(numMolecules, rel_err_tolerance,
                                               autopas::TraversalOption::verletClustersColoring, boxMax);
  }
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, verletClustersColoringTest1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 2e-12;

  for (auto boxMax : {getBoxMaxBig(), getBoxMaxSmall()}) {
    test<autopas::DataLayoutOption::aos, true>(numMolecules, rel_err_tolerance,
                                               autopas::TraversalOption::verletClustersColoring, boxMax);
  }
}

TEST_F(LinkedCellsVersusVerletClusterListsTest, verletClustersColoringTest2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;

  for (auto boxMax : {getBoxMaxBig(), getBoxMaxSmall()}) {
    test<autopas::DataLayoutOption::aos, true>(numMolecules, rel_err_tolerance,
                                               autopas::TraversalOption::verletClustersColoring, boxMax);
  }
}