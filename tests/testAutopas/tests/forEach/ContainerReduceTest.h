/**
 * @file ContainerForEachtest.h
 * @author lgaertner
 * @date 25.08.2021
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/IteratorBehavior.h"
#include "testingHelpers/GenerateValidConfigurations.h"

using testingTuple =
    std::tuple<ContainerConfiguration, bool /*testConstIterators*/, bool /*priorForceCalc*/, autopas::IteratorBehavior>;

class ContainerReduceTest : public AutoPasTestBase, public ::testing::WithParamInterface<testingTuple> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [containerConfig, testConstIterators, priorForceCalc, behavior] = static_cast<ParamType>(info.param);
      std::string str;
      str += containerConfig.container.to_string() + "_";
      str += std::string{"cellSizeFactor"} + std::to_string(containerConfig.cellSizeFactor);
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
   * @tparam AutoPas_T Type of AutoPas container.
   * @param autoPas AutoPas container
   * @return tuple {haloBoxMin, haloBoxMax}
   */
  template <class AutoPas_T>
  auto defaultInit(AutoPas_T &autoPas, const ContainerConfiguration &containerConfig);
};
