/**
 * @file ContainerSwap.h
 * @author F. Gratl
 * @date 14.12.2020
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/LogicHandler.h"
#include "autopas/tuning/selectors/ContainerSelector.h"

namespace autopas {
class ContainerSwapTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<autopas::ContainerOption, autopas::ContainerOption>> {
 public:
  ContainerSwapTest() = default;
  ~ContainerSwapTest() override = default;

  struct twoParamToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      // tuple of ContainerOption
      const auto &[from, to] = static_cast<ParamType>(info.param);
      return "from" + from.to_string() + "To" + to.to_string();
    }
  };

 protected:
  const std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0.1;
  const unsigned int verletRebuildFrequency = 2;
};
}  // namespace autopas
