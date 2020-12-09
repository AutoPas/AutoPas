/**
 * @file FlopCounterTest.h
 * @author F. Gratl
 * @date 01.06.18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

namespace FlopCounterTest {

class FlopCounterTest : public AutoPasTestBase {
 public:
  FlopCounterTest() = default;

  ~FlopCounterTest() override = default;

  void test(autopas::DataLayoutOption dataLayoutOption);
};

}  // end namespace FlopCounterTest
