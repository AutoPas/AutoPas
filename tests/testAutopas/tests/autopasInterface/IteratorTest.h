/**
 * @file IteratorTest.h
 * @author seckler
 * @date 22.07.19
 */

#pragma once

#include <gtest/gtest.h>

#include <tuple>

#include "autopas/AutoPas.h"

namespace IteratorTest {

using testingTuple = std::tuple<autopas::ContainerOption, double /*cell size factor*/, bool /*testConstIterators*/,
                                bool /*priorForceCalc*/>;

class IteratorTest : public testing::Test, public ::testing::WithParamInterface<testingTuple> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = static_cast<ParamType>(info.param);
      std::string str;
      str += containerOption.to_string() + "_";
      str += std::string{"cellSizeFactor"} + std::to_string(cellSizeFactor);
      str += testConstIterators ? "_const" : "_nonConst";
      str += priorForceCalc ? "_priorForceCalc" : "_noPriorForceCalc";
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };

  template <bool testConstIterators>
  static void testOpenMPIterators(autopas::ContainerOption containerOption, double cellSizeFactor,
                                  autopas::IteratorBehavior behavior, bool testRegionIterators, bool priorForceCalc);
};
} // end namespace IteratorTest
