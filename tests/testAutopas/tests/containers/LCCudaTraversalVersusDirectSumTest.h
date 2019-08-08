/**
 * @file LCCudaTraversalVersusDirectSumTest.h
 * @author jspahl
 * @date 11.03.19
 */

#pragma once

#include <gtest/gtest.h>
#include <cstdlib>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "testingHelpers/commonTypedefs.h"

class LCCudaTraversalVersusDirectSumTest : public AutoPasTestBase {
 public:
  LCCudaTraversalVersusDirectSumTest();

  ~LCCudaTraversalVersusDirectSumTest() override = default;

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }

 protected:
  double fRand(double fMin, double fMax) const;

  std::array<double, 3> randomPosition(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) const;

  void fillContainerWithMolecules(
      unsigned long numMolecules,
      autopas::ParticleContainer<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> &cont) const;

  template <bool useNewton3, bool calculateGlobals = false>
  void test(unsigned long numMolecules, double rel_err_tolerance);

  autopas::DirectSum<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> _directSum;
  autopas::LinkedCells<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> _linkedCells;
};
