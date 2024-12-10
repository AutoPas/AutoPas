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
  RebuildNeighborListsTest() = default;
  ~RebuildNeighborListsTest() override = default;

  // Custom function to generate readable names
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[pairwiseConfig, triwiseConfig] = static_cast<ParamType>(info.param);
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
};

inline std::set<autopas::Configuration> getPairwiseConfigs() {
  constexpr double _cellSizeFactor{1.};
  constexpr autopas::Configuration _confVl_list_it_noN3{autopas::ContainerOption::verletLists,
                                                        _cellSizeFactor,
                                                        autopas::TraversalOption::vl_list_iteration,
                                                        autopas::LoadEstimatorOption::none,
                                                        autopas::DataLayoutOption::aos,
                                                        autopas::Newton3Option::disabled,
                                                        autopas::InteractionTypeOption::pairwise};
  constexpr autopas::Configuration _confLc_c01_noN3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c01,        autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise};
  return {_confLc_c01_noN3, _confVl_list_it_noN3};
}

inline std::set<autopas::Configuration> getTriwiseConfigs() {
  constexpr double _cellSizeFactor{1.};
  constexpr autopas::Configuration _confVl_list_it_noN3_3B{autopas::ContainerOption::verletLists,
                                                           _cellSizeFactor,
                                                           autopas::TraversalOption::vl_list_iteration,
                                                           autopas::LoadEstimatorOption::none,
                                                           autopas::DataLayoutOption::aos,
                                                           autopas::Newton3Option::disabled,
                                                           autopas::InteractionTypeOption::triwise};
  constexpr autopas::Configuration _confLc_c01_noN3_3B{
      autopas::ContainerOption::linkedCells,  _cellSizeFactor,
      autopas::TraversalOption::lc_c01,       autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,         autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise};
  constexpr autopas::Configuration _confVl_pairlist_noN3_3B{autopas::ContainerOption::verletLists,
                                                            _cellSizeFactor,
                                                            autopas::TraversalOption::vl_pair_list_iteration_3b,
                                                            autopas::LoadEstimatorOption::none,
                                                            autopas::DataLayoutOption::aos,
                                                            autopas::Newton3Option::disabled,
                                                            autopas::InteractionTypeOption::triwise};
  return {_confLc_c01_noN3_3B, _confVl_list_it_noN3_3B, _confVl_pairlist_noN3_3B};
}
