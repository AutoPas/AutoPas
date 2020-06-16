/**
 * @file TestsAllContainers.h
 * @author humig
 * @date 08.07.2019
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

using ParamType = autopas::ContainerOption;

class AllContainersTests : public AutoPasTestBase, public ::testing::WithParamInterface<ParamType> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      return info.param.to_string();
    };
    return paramToString;
  }

 protected:
  template <class ParticleType = autopas::Particle>
  auto getInitializedContainer() {
    auto containerOptionToTest = GetParam();
    std::array<double, 3> boxMin = {0, 0, 0};
    std::array<double, 3> boxMax = {10, 10, 10};
    double cutoff = 1;
    double skin = 0.2;
    double cellSizeFactor = 1;

    autopas::ContainerSelector<ParticleType, autopas::FullParticleCell<ParticleType>> selector{boxMin, boxMax, cutoff};
    autopas::ContainerSelectorInfo selectorInfo{cellSizeFactor, skin, 32, autopas::LoadEstimatorOption::none};
    selector.selectContainer(containerOptionToTest, selectorInfo);
    return selector.getCurrentContainer();
  }

  void testUpdateContainerDeletesDummy(bool previouslyOwned);
};
