/**
 * @file LinkedCellsVersusVarVerletListsTest.cpp
 * @author humig
 * @date 21.05.19
 *
 * Mostly copied from LinkedCellsVersusVerletLists
 */

#include "LinkedCellsVersusVarVerletListsTest.h"

#include "autopas/containers/verletListsCellBased/verletLists/traversals/VarVerletTraversalAsBuild.h"

LinkedCellsVersusVarVerletListsTest::LinkedCellsVersusVarVerletListsTest()
    : _verletLists(nullptr), _linkedCells(nullptr) {}

template <bool useNewton3, autopas::DataLayoutOption::Value dataLayoutOption>
void LinkedCellsVersusVarVerletListsTest::test(unsigned long numMolecules, double rel_err_tolerance,
                                               std::array<double, 3> boxMax) {
  // generate containers
  _linkedCells = std::make_unique<lctype>(getBoxMin(), boxMax, getCutoff(), 0.1 * getCutoff());
  _verletLists = std::make_unique<vltype>(getBoxMin(), boxMax, getCutoff(), 0.1 * getCutoff(), 4);

  // fill containers
  autopasTools::generators::RandomGenerator::fillWithParticles(
      *_verletLists, Molecule({0., 0., 0.}, {0., 0., 0.}, 0, 0), _verletLists->getBoxMin(), _verletLists->getBoxMax(),
      numMolecules);
  // now fill second container with the molecules from the first one, because
  // otherwise we generate new and different particles
  for (auto it = _verletLists->begin(); it.isValid(); ++it) {
    _linkedCells->addParticle(*it);
  }

  autopas::LJFunctor<Molecule, FMCell> func(getCutoff());
  func.setParticleProperties(1., 1.);
  autopas::VarVerletTraversalAsBuild<FMCell, Molecule, decltype(func), dataLayoutOption, useNewton3> traversalLJV(
      &func);

  autopas::C08Traversal<FMCell, decltype(func), dataLayoutOption, useNewton3> traversalLJ(
      _linkedCells->getCellBlock().getCellsPerDimensionWithHalo(), &func, _linkedCells->getInteractionLength(),
      _linkedCells->getCellBlock().getCellLength());

  _verletLists->rebuildNeighborLists(&traversalLJV);
  _verletLists->iteratePairwise(&traversalLJV);
  _linkedCells->iteratePairwise(&traversalLJ);

  auto itDirect = _verletLists->begin();
  auto itLinked = _linkedCells->begin();

  std::vector<std::array<double, 3>> forcesVerlet(numMolecules), forcesLinked(numMolecules);
  // get and sort by id, the
  for (auto it = _verletLists->begin(); it.isValid(); ++it) {
    Molecule &m = *it;
    forcesVerlet.at(m.getID()) = m.getF();
  }

  for (auto it = _linkedCells->begin(); it.isValid(); ++it) {
    Molecule &m = *it;
    forcesLinked.at(m.getID()) = m.getF();
  }

  for (unsigned long i = 0; i < numMolecules; ++i) {
    for (int d = 0; d < 3; ++d) {
      double f1 = forcesVerlet[i][d];
      double f2 = forcesLinked[i][d];
      EXPECT_NEAR(f1, f2, std::fabs(f1 * rel_err_tolerance));
    }
  }

  autopas::FlopCounterFunctor<Molecule, FMCell> flopsVerlet(getCutoff()), flopsLinked(getCutoff());

  autopas::C08Traversal<FMCell, decltype(flopsLinked), dataLayoutOption, useNewton3> traversalFLOPSLC(
      _linkedCells->getCellBlock().getCellsPerDimensionWithHalo(), &flopsLinked, _linkedCells->getInteractionLength(),
      _linkedCells->getCellBlock().getCellLength());

  autopas::VarVerletTraversalAsBuild<FMCell, Molecule, decltype(flopsLinked), dataLayoutOption, useNewton3>
      traversalFLOPSVerlet(&flopsVerlet);
  _linkedCells->iteratePairwise(&traversalFLOPSLC);
  _verletLists->iteratePairwise(&traversalFLOPSVerlet);

  if (not useNewton3 and dataLayoutOption == autopas::DataLayoutOption::soa) {
    // special case if newton3 is disabled and soa are used: here linked cells will anyways partially use newton3 (for
    // the intra cell interactions), so linked cell kernel calls will be less than for verlet.
    EXPECT_LE(flopsLinked.getKernelCalls(), flopsVerlet.getKernelCalls())
        << "N3: " << (useNewton3 ? "true" : "false") << ", " << autopas::DataLayoutOption(dataLayoutOption).to_string()
        << ", boxMax = [" << boxMax[0] << ", " << boxMax[1] << ", " << boxMax[2] << "]";

  } else {
    // normally the number of kernel calls should be exactly the same
    EXPECT_EQ(flopsLinked.getKernelCalls(), flopsVerlet.getKernelCalls())
        << "N3: " << (useNewton3 ? "true" : "false") << ", " << autopas::DataLayoutOption(dataLayoutOption).to_string()
        << ", boxMax = [" << boxMax[0] << ", " << boxMax[1] << ", " << boxMax[2] << "]";
  }
  // blackbox mode: the following line is only true, if the verlet lists do NOT use less cells than the linked cells
  // (for small scenarios), as the verlet lists fall back to linked cells.
  // Furthermore, if DataLayout=SoA and newton3=disabled, linked cells will still use newton3 for particles inside
  // the same cell. Thus, this case has to be excluded!
  if (dataLayoutOption != autopas::DataLayoutOption::soa and useNewton3) {
    EXPECT_GE(flopsLinked.getDistanceCalculations(), flopsVerlet.getDistanceCalculations());
  }
}

TEST_F(LinkedCellsVersusVarVerletListsTest, test100) {
  unsigned long numMolecules = 100;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-14;

  for (auto boxMax : {std::array<double, 3>{3., 3., 3.}, std::array<double, 3>{10., 10., 10.}}) {
    test<true, autopas::DataLayoutOption::aos>(numMolecules, rel_err_tolerance, boxMax);
    test<true, autopas::DataLayoutOption::soa>(numMolecules, rel_err_tolerance, boxMax);
    test<false, autopas::DataLayoutOption::aos>(numMolecules, rel_err_tolerance, boxMax);
    test<false, autopas::DataLayoutOption::soa>(numMolecules, rel_err_tolerance, boxMax);
  }
}

TEST_F(LinkedCellsVersusVarVerletListsTest, test1000) {
  unsigned long numMolecules = 1000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 2e-12;
  for (auto boxMax : {std::array<double, 3>{3., 3., 3.}, std::array<double, 3>{10., 10., 10.}}) {
    test<true, autopas::DataLayoutOption::aos>(numMolecules, rel_err_tolerance, boxMax);
    test<true, autopas::DataLayoutOption::soa>(numMolecules, rel_err_tolerance, boxMax);
    test<false, autopas::DataLayoutOption::aos>(numMolecules, rel_err_tolerance, boxMax);
    test<false, autopas::DataLayoutOption::soa>(numMolecules, rel_err_tolerance, boxMax);
  }
}

TEST_F(LinkedCellsVersusVarVerletListsTest, test2000) {
  unsigned long numMolecules = 2000;

  // empirically determined and set near the minimal possible value
  // i.e. if something changes, it may be needed to increase value
  // (and OK to do so)
  double rel_err_tolerance = 1e-10;
  for (auto boxMax : {std::array<double, 3>{3., 3., 3.}, std::array<double, 3>{10., 10., 10.}}) {
    test<true, autopas::DataLayoutOption::aos>(numMolecules, rel_err_tolerance, boxMax);
    test<true, autopas::DataLayoutOption::soa>(numMolecules, rel_err_tolerance, boxMax);
    test<false, autopas::DataLayoutOption::aos>(numMolecules, rel_err_tolerance, boxMax);
    test<false, autopas::DataLayoutOption::soa>(numMolecules, rel_err_tolerance, boxMax);
  }
}
