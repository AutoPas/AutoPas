/**
 * @file RebuildNeighborListsTest.h
 * @author muehlhaeusser
 * @date 29.11.2024
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/tuning/Configuration.h"

class RebuildNeighborListsTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<autopas::Configuration, autopas::Configuration>> {
 public:
  static std::set<autopas::Configuration> getPairwiseConfigs();
  static std::set<autopas::Configuration> getTriwiseConfigs();

  // Custom function to generate readable names
  static std::string configsToString(
      const ::testing::TestParamInfo<std::tuple<autopas::Configuration, autopas::Configuration>> &info) {
    const auto &configs = info.param;
    const auto &pairwiseConfig = std::get<0>(info.param);
    const auto &triwiseConfig = std::get<1>(info.param);
    std::stringstream resStream;
    resStream << "2B_" << pairwiseConfig.container.to_string() << "_" << pairwiseConfig.traversal.to_string() << "_"
              << pairwiseConfig.dataLayout.to_string() << "_"
              << (pairwiseConfig.newton3 == autopas::Newton3Option::enabled ? "_N3" : "_noN3");
    resStream << "__3B_" << triwiseConfig.container.to_string() << "_" << triwiseConfig.traversal.to_string() << "_"
              << triwiseConfig.dataLayout.to_string() << "_"
              << (triwiseConfig.newton3 == autopas::Newton3Option::enabled ? "_N3" : "_noN3");
    std::string res = resStream.str();
    std::replace(res.begin(), res.end(), '-', '_');
    std::replace(res.begin(), res.end(), '.', '_');
    return res;
  }
};
