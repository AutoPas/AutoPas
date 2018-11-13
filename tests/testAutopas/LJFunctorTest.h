/**
 * @file LJFunctorTest.h
 * @author seckler
 * @date 06.11.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "testingHelpers/RandomGenerator.h"

class LJFunctorTest : public AutoPasTestBase {
 public:
  LJFunctorTest() : AutoPasTestBase() {
    cutoff = 1.;
    epsilon = 1.;
    sigma = 1.;
    shift = 0.1;
    lowCorner = {0., 0., 0.};
    highCorner = {5., 5., 5.};
    expectedForce = {-4547248.8989645941, -9094497.7979291882, -13641746.696893783};
    expectedEnergy = 3178701.6514326506 / 6.;
    expectedVirial = 6366148.4585504318;
    absDelta = 1e-7;
  }

  void SetUp() override{};

  void TearDown() override{};

 protected:
  void testAoSNoGlobals(bool newton3);

  enum where_type { inside, boundary, outside };
  void testAoSGlobals(where_type where, bool newton3, bool duplicatedCalculation);

  double cutoff;
  double epsilon;
  double sigma;
  double shift;
  std::array<double, 3> lowCorner;
  std::array<double, 3> highCorner;

  std::array<double, 3> expectedForce;

  double expectedVirial;
  double expectedEnergy;
  double absDelta;
};
