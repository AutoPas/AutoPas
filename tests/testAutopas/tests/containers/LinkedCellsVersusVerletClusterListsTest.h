/**
 * @file LinkedCellsVersusVerletClusterListsTest.h
 * @author nguyen
 * @date 21.10.18
 */

#pragma once

#include <gtest/gtest.h>
#include <cstdlib>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/molecularDynamics/ParticleClassLibrary.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsVersusVerletClusterListsTest : public AutoPasTestBase {
 public:
  ~LinkedCellsVersusVerletClusterListsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMaxSmall() const { return {3.0, 3.0, 3.0}; }
  std::array<double, 3> getBoxMaxBig() const { return {10.0, 10.0, 10.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  template <autopas::DataLayoutOption dataLayout, bool useNewton3>
  void test(unsigned long numMolecules, double rel_err_tolerance, autopas::TraversalOption traversalOption,
            std::array<double, 3> boxMax);

  using Verlet = autopas::VerletClusterLists<autopas::MoleculeLJ<>>;
  using Linked = autopas::LinkedCells<autopas::MoleculeLJ<>, autopas::FullParticleCell<autopas::MoleculeLJ<>>>;
};
