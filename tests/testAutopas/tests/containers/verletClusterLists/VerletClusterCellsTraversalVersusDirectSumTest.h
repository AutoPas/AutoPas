/**
 * @file VerletClusterCellsTraversalVersusDirectSumTest.h
 * @author jspahl
 * @date 3.04.19
 */

#pragma once

#include <gtest/gtest.h>
#include <cstdlib>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/containers/verletClusterLists/VerletClusterCells.h"
#include "testingHelpers/commonTypedefs.h"

class VerletClusterCellsTraversalVersusDirectSumTest : public AutoPasTestBase {
 public:
  VerletClusterCellsTraversalVersusDirectSumTest();

  ~VerletClusterCellsTraversalVersusDirectSumTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  double fRand(double fMin, double fMax) const;

  std::array<double, 3> randomPosition(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) const;

  void fillContainerWithMolecules(unsigned long numMolecules, autopas::ParticleContainer<FMCell> &cont) const;

  template <bool useNewton3, autopas::DataLayoutOption dataLayout = autopas::DataLayoutOption::aos,
            bool calculateGlobals = false>
  void test(unsigned long numMolecules, double rel_err_tolerance);

  autopas::DirectSum<FMCell> _directSum;
  autopas::VerletClusterCells<Molecule> _verletCluster;
};
