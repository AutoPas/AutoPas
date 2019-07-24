/**
 * @file IteratorTest.h
 * @author seckler
 * @date 22.07.19
 */

#pragma once

#include <gtest/gtest.h>
#include <tuple>
#include "autopas/AutoPas.h"

using testingTuple = std::tuple<autopas::ContainerOption, double /*cell size factor*/>;

class IteratorTest : public testing::Test, public ::testing::WithParamInterface<testingTuple> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string str;
      str += autopas::utils::StringUtils::to_string(std::get<0>(inputTuple)) + "_";
      str += std::string{"cellSizeFactor"} + autopas::utils::StringUtils::to_string(std::get<1>(inputTuple));
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };
};