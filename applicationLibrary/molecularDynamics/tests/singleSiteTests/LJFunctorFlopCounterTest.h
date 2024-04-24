/**
 * @file LJFunctorFlopCounterTest.h
 * @author F. Gratl
 * @date 01.06.18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/DataLayoutOption.h"

class LJFunctorFlopCounterTest : public AutoPasTestBase {
 public:
  LJFunctorFlopCounterTest() = default;

  ~LJFunctorFlopCounterTest() override = default;

  void testFLOPCounter(autopas::DataLayoutOption dataLayoutOption, bool newton3);
};
