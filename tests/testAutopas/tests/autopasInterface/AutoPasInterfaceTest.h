/**
 * @file AutoPasInterfaceTest.h
 * @author seckler
 * @date 13.05.19
 */

#pragma once

#include <gtest/gtest.h>
#include "autopas/AutoPas.h"

class AutoPasInterfaceTest : public testing::Test, public ::testing::WithParamInterface<autopas::ContainerOption> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string containerSelector_str(autopas::utils::StringUtils::to_string(inputTuple));
      return containerSelector_str;
    }
  };
};

class AutoPasInterface2ContainersTest
    : public testing::Test,
      public ::testing::WithParamInterface<std::tuple<autopas::ContainerOption, autopas::ContainerOption>> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string from_str(autopas::utils::StringUtils::to_string(std::get<0>(inputTuple)));
      std::string to_str(autopas::utils::StringUtils::to_string(std::get<1>(inputTuple)));
      // replace all '-' with '_', otherwise the test name is invalid
      // std::replace(traversal.begin(), traversal.end(), '-', '_');
      return from_str + "And" + to_str;
    }
  };
};
