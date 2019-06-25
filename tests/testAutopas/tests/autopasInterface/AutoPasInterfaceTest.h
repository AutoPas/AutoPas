/**
 * @file AutoPasInterfaceTest.h
 * @author seckler
 * @date 13.05.19
 */

#pragma once

#include <gtest/gtest.h>
#include <tuple>
#include "autopas/AutoPas.h"

using testingTuple = std::tuple<std::tuple<autopas::ContainerOption, autopas::TraversalOption>,
                                autopas::DataLayoutOption, autopas::Newton3Option, double /*cell size factor*/>;

class AutoPasInterfaceTest : public testing::Test, public ::testing::WithParamInterface<testingTuple> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string str;
      str += autopas::utils::StringUtils::to_string(std::get<0>(std::get<0>(inputTuple))) + "_";
      str += autopas::utils::StringUtils::to_string(std::get<1>(std::get<0>(inputTuple))) + "_";
      str += autopas::utils::StringUtils::to_string(std::get<1>(inputTuple)) + "_";
      str += autopas::utils::StringUtils::to_string(std::get<2>(inputTuple)) + "_";
      str += std::string{"cellSizeFactor"} + autopas::utils::StringUtils::to_string(std::get<3>(inputTuple));
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
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
