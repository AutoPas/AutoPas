/**
 * @file AutoPasInterfaceTest.h
 * @author seckler
 * @date 13.05.19
 */

#pragma once

#include <gtest/gtest.h>

#include <tuple>

#include "autopas/AutoPas.h"

namespace AutoPasInterfaceTest {

using testingTuple =
    std::tuple<std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::LoadEstimatorOption>,
               autopas::DataLayoutOption, autopas::Newton3Option, double /*cell size factor*/>;

class AutoPasInterfaceTest : public testing::Test, public ::testing::WithParamInterface<testingTuple> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string str;
      str += std::get<0>(std::get<0>(inputTuple)).to_string() + "_";
      str += std::get<1>(std::get<0>(inputTuple)).to_string() + "_";
      str += std::get<2>(std::get<0>(inputTuple)).to_string() + "_";
      str += std::get<1>(inputTuple).to_string() + "_";
      str += "N3" + std::get<2>(inputTuple).to_string() + "_";
      str += std::string{"cellSizeFactor"} + std::to_string(std::get<3>(inputTuple));
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
      std::string fromStr(std::get<0>(inputTuple).to_string());
      std::string toStr(std::get<1>(inputTuple).to_string());
      // replace all '-' with '_', otherwise the test name is invalid
      // std::replace(traversal.begin(), traversal.end(), '-', '_');
      return fromStr + "And" + toStr;
    }
  };
};

}  // end namespace AutoPasInterfaceTest
