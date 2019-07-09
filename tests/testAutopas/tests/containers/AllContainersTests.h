/**
 * @file TestsAllContainers.h
 * @author humig
 * @date 08.07.2019
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

using ParamType = autopas::ContainerOption;

class AllContainersTests : public AutoPasTestBase, public ::testing::WithParamInterface<ParamType> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      return autopas::utils::StringUtils::to_string(info.param);
    };
    return paramToString;
  }

 protected:
  void SetUp() override {
    auto containerOptionToTest = GetParam();
    std::array<double, 3> boxMin = {0, 0, 0};
    std::array<double, 3> boxMax = {10, 10, 10};
    double cutoff = 1;
    double skin = 0.2;
    double cellSizeFactor = 1;

    autopas::ContainerSelector<Particle, FPCell> selector{boxMin, boxMax, cutoff};
    autopas::ContainerSelectorInfo selectorInfo{cellSizeFactor, skin};
    selector.selectContainer(containerOptionToTest, selectorInfo);
    _container = selector.getCurrentContainer();
  }

 protected:
  std::shared_ptr<autopas::ParticleContainer<Particle, FPCell>> _container{};
};
