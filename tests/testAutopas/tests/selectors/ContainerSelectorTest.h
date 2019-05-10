/**
 * @file ContainerSelectorTest.h
 * @author F. Gratl
 * @date 22.06.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

class ContainerSelectorTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<autopas::ContainerOption, autopas::ContainerOption>> {
 public:
  ContainerSelectorTest() = default;
  ~ContainerSelectorTest() override = default;

  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string from_str(autopas::utils::StringUtils::to_string(std::get<0>(inputTuple)));
      std::string to_str(autopas::utils::StringUtils::to_string(std::get<1>(inputTuple)));
      // replace all '-' with '_', otherwise the test name is invalid
      // std::replace(traversal.begin(), traversal.end(), '-', '_');
      return "from" + from_str + "To" + to_str;
    }
  };
};
