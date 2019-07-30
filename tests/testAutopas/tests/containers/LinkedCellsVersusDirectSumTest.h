/**
 * @file LinkedCellsVersusDirectSumTest.h
 * @author tchipev
 * @date 23.01.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/molecularDynamics/ParticleClassLibrary.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsVersusDirectSumTest : public AutoPasTestBase {
 public:
  LinkedCellsVersusDirectSumTest();

  ~LinkedCellsVersusDirectSumTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  void test(unsigned long numMolecules, double rel_err_tolerance);

  autopas::DirectSum<autopas::MoleculeLJ<>, autopas::FullParticleCell<autopas::MoleculeLJ<>>> _directSum;
  autopas::LinkedCells<autopas::MoleculeLJ<>, autopas::FullParticleCell<autopas::MoleculeLJ<>>> _linkedCells;
};
