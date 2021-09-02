/**
 * @file ContainerForEachtest.h
 * @author lgaertner
 * @date 25.08.2021
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/IteratorBehavior.h"

using testingTuple = std::tuple<autopas::ContainerOption, double /*cell size factor*/, bool /*testConstIterators*/,
    bool /*priorForceCalc*/, autopas::IteratorBehavior>;

class ContainerReduceTest : public AutoPasTestBase, public ::testing::WithParamInterface<testingTuple> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc, behavior] =
      static_cast<ParamType>(info.param);
      std::string str;
      str += containerOption.to_string() + "_";
      str += std::string{"cellSizeFactor"} + std::to_string(cellSizeFactor);
      str += testConstIterators ? "_const" : "_nonConst";
      str += priorForceCalc ? "_priorForceCalc" : "_noPriorForceCalc";
      str += "_" + behavior.to_string();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };

  /**
   * Initialize the given AutoPas object with the default values for this test class.
   * @tparam AutoPasT
   * @param autoPas
   * @return tuple {haloBoxMin, haloBoxMax}
   */
  template <class AutoPasT>
  auto defaultInit(AutoPasT &autoPas, autopas::ContainerOption &containerOption, double cellSizeFactor);
};
