/**
 * @file AllContainersTests.h
 * @author humig
 * @date 08.07.2019
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

class AllContainersTestsBase : public AutoPasTestBase {
 protected:
  std::array<double, 3> boxMin = {0, 0, 0};
  std::array<double, 3> boxMax = {10, 10, 10};

  template <class ParticleType = autopas::Particle>
  auto getInitializedContainer(autopas::ContainerOption containerOptionToTest) {
    double cutoff = 1;
    double skin = 0.2;
    double cellSizeFactor = 1;

    autopas::ContainerSelector<Particle> selector{boxMin, boxMax, cutoff};
    autopas::ContainerSelectorInfo selectorInfo{cellSizeFactor, skin, 32, autopas::LoadEstimatorOption::none};
    selector.selectContainer(containerOptionToTest, selectorInfo);
    return selector.getCurrentContainer();
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