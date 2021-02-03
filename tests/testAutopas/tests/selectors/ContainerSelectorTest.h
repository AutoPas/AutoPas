/**
 * @file ContainerSelectorTest.h
 * @author F. Gratl
 * @date 22.06.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

class ContainerSelectorTest : public AutoPasTestBase, public ::testing::WithParamInterface<autopas::ContainerOption> {
 public:
  ContainerSelectorTest() = default;
  ~ContainerSelectorTest() override = default;

  struct oneParamToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      // tuple of ContainerOption
      const auto &option = static_cast<ParamType>(info.param);
      return option.to_string();
    }
  };

 protected:
  const std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0.1;
};
