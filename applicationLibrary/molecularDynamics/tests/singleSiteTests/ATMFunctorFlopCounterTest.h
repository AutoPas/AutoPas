/**
 * @file ATMFunctorFlopCounterTest.h
 * @author muehlhaeusser
 * @date 30.07.24
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/DataLayoutOption.h"

class ATMFunctorFlopCounterTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<autopas::DataLayoutOption, bool, bool, bool>> {
 public:
  ATMFunctorFlopCounterTest() = default;

  ~ATMFunctorFlopCounterTest() override = default;

  template <bool calculateGlobals>
  void testFLOPCounter(autopas::DataLayoutOption dataLayoutOption, bool newton3, bool isVerlet);

  template <bool calculateGlobals>
  void testFLOPCounterAoSOMP(bool newton3);

  template <bool calculateGlobals>
  void testFLOPCounterSoASingleAndPairOMP(bool newton3){};

  template <bool calculateGlobals>
  void testFLOPCounterSoAVerletOMP(bool newton3){};
};
