/**
 * @file LinkedCellsVersusDirectSumTest.h
 * @author seckler
 * @date 21.05.18
 */

#pragma once

#include <gtest/gtest.h>

#include <cstdlib>
#include <memory>

#include "AutoPasTestBase.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsVersusVerletListsTest : public AutoPasTestBase {
 public:
  LinkedCellsVersusVerletListsTest();

  ~LinkedCellsVersusVerletListsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  double getCutoff() const { return .9; }

 protected:
  template <bool useNewton3, autopas::DataLayoutOption::Value dataLayoutOption>
  void test(unsigned long numMolecules, double rel_err_tolerance, std::array<double, 3> boxMax);

  using vltype = autopas::VerletLists<Molecule>;
  using lctype = autopas::LinkedCells<FMCell>;
  std::unique_ptr<vltype> _verletLists;
  std::unique_ptr<lctype> _linkedCells;
};
