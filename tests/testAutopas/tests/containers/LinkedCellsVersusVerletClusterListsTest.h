/**
 * @file LinkedCellsVersusVerletClusterListsTest.h
 * @author nguyen
 * @date 21.10.18
 */

#pragma once

#include <gtest/gtest.h>
#include <cstdlib>
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class LinkedCellsVersusVerletClusterListsTest : public AutoPasTestBase {
 public:
  ~LinkedCellsVersusVerletClusterListsTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  template <autopas::DataLayoutOption dataLayout, bool useNewton3>
  void test(unsigned long numMolecules, double rel_err_tolerance);

  using Verlet = autopas::VerletClusterLists<autopas::MoleculeLJ>;
  using Linked = autopas::LinkedCells<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>;
};
