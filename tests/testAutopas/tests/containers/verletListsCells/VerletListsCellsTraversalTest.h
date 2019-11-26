/**
 * @file VerletListsCellsTraversalTest.h
 * @author nguyen
 * @date 26.09.18
 */

#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstdlib>

#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class VerletListsCellsTraversalTest : public AutoPasTestBase {
 public:
  VerletListsCellsTraversalTest();

  ~VerletListsCellsTraversalTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  void test(unsigned long numMolecules);

  autopas::VerletListsCells<Molecule> _verletListsCells, _verletListsCells_cs2;
};
