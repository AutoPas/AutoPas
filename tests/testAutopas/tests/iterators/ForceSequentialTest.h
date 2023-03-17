/**
 * @file ForceSequentialTest.h
 * @author F. Gratl
 * @date 15.02.2023
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"

using testingTuple = std::tuple<autopas::ContainerOption>;

class ForceSequentialTest : public AutoPasTestBase, public ::testing::WithParamInterface<testingTuple> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc, behavior] =
          static_cast<ParamType>(info.param);
      std::string str;
      str += containerOption.to_string();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };
};
