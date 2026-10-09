/**
 * @file AllContainersTests.h
 * @author humig
 * @date 08.07.2019
 */

#pragma once

#include <algorithm>
#include <string>

#include "AutoPasTestBase.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "testingHelpers/GenerateValidConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

class AllContainersTestsBase : public AutoPasTestBase {
 protected:
  std::array<double, 3> boxMin = {0, 0, 0};
  std::array<double, 3> boxMax = {10, 10, 10};
  double cutoff = 1;
  const double skin = 0.2;
  const unsigned int rebuildFrequency = 20;

  template <class Particle_T>
  auto getInitializedContainer(const ContainerConfiguration &containerConfig) {
    const autopas::ContainerSelectorInfo selectorInfo{
        boxMin, boxMax, cutoff, containerConfig.cellSizeFactor, skin, 32, 8, 8, autopas::LoadEstimatorOption::none};
    auto container =
        autopas::ContainerSelector<Particle_T>::generateContainer(containerConfig.container, selectorInfo);
    return std::move(container);
  }
};

using ParamType = ContainerConfiguration;

class AllContainersTests : public AllContainersTestsBase, public ::testing::WithParamInterface<ParamType> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      const auto &containerConfig = info.param;
      return containerConfig.container.to_string() + "_cellSizeFactor" +
             std::to_string(containerConfig.cellSizeFactor);
    };
    return paramToString;
  }
};

using ParamTypeBothUpdates = std::tuple<ContainerConfiguration, bool /*keep Lists Valid*/>;

class AllContainersTestsBothUpdates : public AllContainersTestsBase,
                                      public ::testing::WithParamInterface<ParamTypeBothUpdates> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      auto [containerConfig, keepListValid] = info.param;
      std::string str = containerConfig.container.to_string() + "_cellSizeFactor" +
                        std::to_string(containerConfig.cellSizeFactor) + "_" +
                        (keepListValid ? "keepListsValid" : "allowListInvalidation");
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    };
    return paramToString;
  }

 protected:
  void testUpdateContainerDeletesDummy(bool previouslyOwned);
};