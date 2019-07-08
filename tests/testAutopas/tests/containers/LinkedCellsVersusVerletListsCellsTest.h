/**
 * @file LinkedCellsVersusVerletListsCellsTest.h
 * @author nguyen
 * @date 08.09.18
 */

#pragma once

#include <gtest/gtest.h>
#include <cstdlib>
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsVersusVerletListsCellsTest : public AutoPasTestBase {
 public:
  LinkedCellsVersusVerletListsCellsTest();

  ~LinkedCellsVersusVerletListsCellsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  void test(unsigned long numMolecules, double rel_err_tolerance);

  autopas::VerletListsCells<autopas::MoleculeLJ> _verletListsCells;
  autopas::LinkedCells<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> _linkedCells;
};
