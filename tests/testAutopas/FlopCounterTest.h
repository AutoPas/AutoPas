/**
 * FlopCounterTest.h
 *
 *  Created on: 6/1/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <testingHelpers/commonTypedefs.h>
#include "AutoPas.h"
#include "AutoPasTestBase.h"

class FlopCounterTest : public AutoPasTestBase {
 public:
  FlopCounterTest() = default;

  ~FlopCounterTest() override = default;

  void test(autopas::DataLayoutOption dataLayoutOption);
};
