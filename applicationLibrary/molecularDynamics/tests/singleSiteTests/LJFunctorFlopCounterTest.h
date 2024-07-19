/**
 * @file LJFunctorFlopCounterTest.h
 * @author F. Gratl
 * @date 01.06.18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/DataLayoutOption.h"

class LJFunctorFlopCounterTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<autopas::DataLayoutOption, bool, bool, bool, bool, bool>> {
 public:
  LJFunctorFlopCounterTest() = default;

  ~LJFunctorFlopCounterTest() override = default;

  template <bool useMixing, bool calculateGlobals, bool applyShift>
  void testFLOPCounter(autopas::DataLayoutOption dataLayoutOption, bool newton3, bool isVerlet);

  template <bool useMixing, bool calculateGlobals, bool applyShift>
  void testFLOPCounterAoSOMP(bool newton3);

  template <bool useMixing, bool calculateGlobals, bool applyShift>
  void testFLOPCounterSoASingleAndPairOMP(bool newton3);

  template <bool useMixing, bool calculateGlobals, bool applyShift>
  void testFLOPCounterSoAVerletOMP(bool newton3);
};
