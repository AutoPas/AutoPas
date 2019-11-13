/**
 * @file LinkedCellsVersusVarVerletListsTest.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include <cstdlib>
#include <memory>

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/containers/verletListsCellBased/verletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "testingHelpers/RandomGenerator.h"
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

  using vltype = autopas::VarVerletLists<autopas::MoleculeLJ, autopas::VerletNeighborListAsBuild<autopas::MoleculeLJ>>;
  using lctype = autopas::LinkedCells<autopas::FullParticleCell<autopas::MoleculeLJ>>;
  std::unique_ptr<vltype> _verletLists;
  std::unique_ptr<lctype> _linkedCells;
};
