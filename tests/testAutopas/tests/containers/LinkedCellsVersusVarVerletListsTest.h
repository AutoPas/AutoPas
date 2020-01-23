/**
 * @file LinkedCellsVersusVarVerletListsTest.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include <gtest/gtest.h>

#include <cstdlib>
#include <memory>

#include "AutoPasTestBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsVersusVarVerletListsTest : public AutoPasTestBase {
 public:
  LinkedCellsVersusVarVerletListsTest();

  ~LinkedCellsVersusVarVerletListsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  double getCutoff() const { return .9; }

 protected:
  template <bool useNewton3, autopas::DataLayoutOption::Value dataLayoutOption>
  void test(unsigned long numMolecules, double rel_err_tolerance, std::array<double, 3> boxMax);

  using vltype = autopas::VarVerletLists<Molecule, autopas::VerletNeighborListAsBuild<Molecule>>;
  using lctype = autopas::LinkedCells<FMCell>;
  std::unique_ptr<vltype> _verletLists;
  std::unique_ptr<lctype> _linkedCells;
};
