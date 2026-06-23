/**
 * @file RegionParticleIteratorTest.h
 * @author seckler
 * @date 03.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/IteratorBehavior.h"
#include "testingHelpers/GenerateValidConfigurations.h"

using testingTupleOne = std::tuple<ContainerConfiguration, bool /*testConstIterators*/,
                                   bool /*priorForceCalc*/, autopas::IteratorBehavior>;

class RegionParticleIteratorTestBase : public AutoPasTestBase {
 protected:
  /**
   * Initialize the given AutoPas object with the default values for this test class.
   * @tparam AutoPasT
   * @param autoPas
   * @param containerConfig container-CSF pair
   * @return tuple {haloBoxMin, haloBoxMax}
   */
  template <class AutoPasT>
  auto defaultInit(AutoPasT &autoPas, ContainerConfiguration containerConfig);
};

class RegionParticleIteratorTestOne : public RegionParticleIteratorTestBase,
                                      public ::testing::WithParamInterface<testingTupleOne> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [containerConfig, testConstIterators, priorForceCalc, behavior] = static_cast<ParamType>(info.param);
      std::string str;
      str += containerConfig.container.to_string() + "_";
      str += std::string{"csf"} + "_" + std::to_string(containerConfig.cellSizeFactor);
      str += testConstIterators ? "_const" : "_nonConst";
      str += priorForceCalc ? "_priorForceCalc" : "_noPriorForceCalc";
      str += "_" + behavior.to_string();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };
};

class RegionParticleIteratorTestTwo : public RegionParticleIteratorTestBase,
                                      public ::testing::WithParamInterface<ContainerConfiguration> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto containerOption = static_cast<ParamType>(info.param);
      std::string str;
      str += containerOption.container.to_string() + "_";
      str += std::string{"csf"} + "_" + std::to_string(containerOption.cellSizeFactor);
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };
};