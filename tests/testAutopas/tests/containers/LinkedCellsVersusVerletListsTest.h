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
#include "autopas/autopasIncludes.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsVersusVerletListsTest : public AutoPasTestBase {
 public:
  LinkedCellsVersusVerletListsTest();

  ~LinkedCellsVersusVerletListsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  double getCutoff() const { return .9; }

 protected:
  void test(unsigned long numMolecules, double rel_err_tolerance, std::array<double, 3> boxMax, bool useSoA,
            bool blackBoxMode);

  using vltype = autopas::VerletLists<autopas::MoleculeLJ>;
  using lctype = autopas::LinkedCells<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>;
  std::unique_ptr<vltype> _verletLists;
  std::unique_ptr<lctype> _linkedCells;
};
