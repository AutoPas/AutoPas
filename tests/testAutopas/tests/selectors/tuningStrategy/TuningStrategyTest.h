/**
 * @file TuningStrategyTest.h
 * @author Jan Nguyen
 * @date 27.04.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/selectors/tuningStrategy/TuningStrategyFactory.h"

namespace TuningStrategyTest {

class TuningStrategyTest : public AutoPasTestBase, public ::testing::WithParamInterface<autopas::TuningStrategyOption> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto tuningStrategy(info.param.to_string());
      // replace all '-' with '_', otherwise the test name is invalid
      std::replace(tuningStrategy.begin(), tuningStrategy.end(), '-', '_');
      return tuningStrategy;
    }
  };
};

} // end namespace TuningStrategyTest
