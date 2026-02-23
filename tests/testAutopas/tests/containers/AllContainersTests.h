/**
 * @file AllContainersTests.h
 * @author humig
 * @date 08.07.2019
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

class AllContainersTestsBase : public AutoPasTestBase {
 protected:
  std::array<double, 3> boxMin = {0, 0, 0};
  std::array<double, 3> boxMax = {10, 10, 10};
  double cutoff = 1;
  const double skin = 0.2;
  const unsigned int rebuildFrequency = 20;
  const double cellSizeFactor = 1;

  template <class Particle_T>
  auto getInitializedContainer(autopas::ContainerOption containerOptionToTest) {
    const autopas::ContainerSelectorInfo selectorInfo{boxMin, boxMax, cutoff, cellSizeFactor,
                                                      skin,   32,     8,      autopas::LoadEstimatorOption::none};
    auto container = autopas::ContainerSelector<Particle_T>::generateContainer(containerOptionToTest, selectorInfo);
    return std::move(container);
  }
};

using ParamType = std::tuple<autopas::ContainerOption>;

class AllContainersTests : public AllContainersTestsBase, public ::testing::WithParamInterface<ParamType> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      auto [containerOption] = info.param;
      return containerOption.to_string();
    };
    return paramToString;
  }
};

using ParamTypeBothUpdates = std::tuple<autopas::ContainerOption, bool /*keep Lists Valid*/>;

class AllContainersTestsBothUpdates : public AllContainersTestsBase,
                                      public ::testing::WithParamInterface<ParamTypeBothUpdates> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      auto [containerOption, keepListValid] = info.param;
      return containerOption.to_string() + "_" + (keepListValid ? "keepListsValid" : "allowListInvalidation");
    };
    return paramToString;
  }

 protected:
  void testUpdateContainerDeletesDummy(bool previouslyOwned);
};