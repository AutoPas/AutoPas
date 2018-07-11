/**
 * @file LinkedCellsVersusDirectSumTest.h
 * @author seckler
 * @date 21.05.18
 */

#pragma once

#include <gtest/gtest.h>
#include <cstdlib>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "testingHelpers/RandomGenerator.h"

class LinkedCellsVersusVerletListsTest : public AutoPasTestBase {
 public:
  LinkedCellsVersusVerletListsTest();

  ~LinkedCellsVersusVerletListsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  void test(unsigned long numMolecules, double rel_err_tolerance);

  autopas::VerletLists<autopas::MoleculeLJ> _verletLists;
  autopas::LinkedCells<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> _linkedCells;
};
