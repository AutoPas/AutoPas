/**
 * @file FlopCounterTest.h
 * @author F. Gratl
 * @date 01.06.18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/DataLayoutOption.h"

class FlopCounterTest : public AutoPasTestBase {
 public:
  FlopCounterTest() = default;

  ~FlopCounterTest() override = default;

  void test(autopas::DataLayoutOption dataLayoutOption);
};
