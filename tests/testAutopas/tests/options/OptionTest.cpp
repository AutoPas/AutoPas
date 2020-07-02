/**
 * @file OptionTest.cpp
 * @author F. Gratl
 * @date 29.10.2019
 */

#include "OptionTest.h"

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "tests/utils/StringUtilsTest.h"

// parseOptions tests
// these tests shall not (yet :) be generated since we here want to pass strings that do not match exactly.

/**
 * The parseXXXOptionsTests define a mapping of enums to strings. It is then tested if parseOptions
 * can correctly parse them individually or all at once.
 */

TEST(OptionTest, parseTraversalOptionsTest) {
  std::map<autopas::TraversalOption, std::string> mapEnumString = {
      {autopas::TraversalOption::lc_c01, "lc_c01"},
      {autopas::TraversalOption::vlc_c01, "verlec01"},
      {autopas::TraversalOption::lc_c01_combined_SoA, "c01-combined"},
      {autopas::TraversalOption::lc_c01_cuda, "cudac01"},
      {autopas::TraversalOption::lc_c04, "lc_c04"},
      {autopas::TraversalOption::lc_c04_combined_SoA, "c04-soa"},
      {autopas::TraversalOption::lc_c08, "lc_c08"},
      {autopas::TraversalOption::lc_c18, "lc_c18"},
      {autopas::TraversalOption::vlc_c18, "verletc18"},
      {autopas::TraversalOption::ds_sequential, "direct"},
      {autopas::TraversalOption::lc_sliced, "slicedv01"},
      {autopas::TraversalOption::vlc_sliced, "verlet-sliced"},
      {autopas::TraversalOption::vvl_as_built, "var-verlet-lists-as-build"},
      {autopas::TraversalOption::vcc_cluster_iteration, "verlet-cluster-cells"},
      {autopas::TraversalOption::vcl_cluster_iteration, "verlet-clusters"},
      {autopas::TraversalOption::vcl_c06, "verlet-clusters-coloring"},
      {autopas::TraversalOption::vl_list_iteration, "verletlists"},
      {autopas::TraversalOption::vcl_c01_balanced, "verlet-cluste-static"},
      {autopas::TraversalOption::lc_sliced_balanced, "balancedsliced"},
      {autopas::TraversalOption::vlc_sliced_balanced, "balancedsliced-verlet"},
      {autopas::TraversalOption::lc_c04_HCP, "lc_c04_HCP"},
  };

  EXPECT_EQ(mapEnumString.size(), autopas::TraversalOption::getOptionNames().size());

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

TEST(OptionTest, parseContainerOptionsTest) {
  std::map<autopas::ContainerOption, std::string> mapEnumString = {
      {autopas::ContainerOption::directSum, "directSum"},
      {autopas::ContainerOption::linkedCells, "linkedCells"},
      {autopas::ContainerOption::varVerletListsAsBuild, "varVerletListsAsBuild"},
      {autopas::ContainerOption::verletClusterCells, "vclustercells"},
      {autopas::ContainerOption::verletClusterLists, "vclusterlists"},
      {autopas::ContainerOption::verletLists, "verletLists"},
      {autopas::ContainerOption::verletListsCells, "verletLists-cells"},
  };

  EXPECT_EQ(mapEnumString.size(), autopas::ContainerOption::getOptionNames().size());

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

TEST(OptionTest, parseDataLayoutOptionsTest) {
  std::map<autopas::DataLayoutOption, std::string> mapEnumString = {
    {autopas::DataLayoutOption::aos, "aos"},
    {autopas::DataLayoutOption::soa, "soa"},
#if defined(AUTOPAS_CUDA)
    {autopas::DataLayoutOption::cuda, "cuda"},
#endif
  };

  EXPECT_EQ(mapEnumString.size(), autopas::DataLayoutOption::getOptionNames().size());

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

TEST(OptionTest, parseSelectorOptionsTest) {
  std::map<autopas::SelectorStrategyOption, std::string> mapEnumString = {
      {autopas::SelectorStrategyOption::fastestAbs, "absolute"},
      {autopas::SelectorStrategyOption::fastestMean, "mean"},
      {autopas::SelectorStrategyOption::fastestMedian, "median"},
  };

  EXPECT_EQ(mapEnumString.size(), autopas::SelectorStrategyOption::getOptionNames().size());

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

TEST(OptionTest, parseTuningStrategyOptionsTest) {
  std::map<autopas::TuningStrategyOption, std::string> mapEnumString = {
      {autopas::TuningStrategyOption::bayesianSearch, "bayesian"},
      {autopas::TuningStrategyOption::bayesianClusterSearch, "bayesian-cluster"},
      {autopas::TuningStrategyOption::fullSearch, "full"},
      {autopas::TuningStrategyOption::randomSearch, "random"},
      {autopas::TuningStrategyOption::activeHarmony, "harmony"},
      {autopas::TuningStrategyOption::predictiveTuning, "predictive"},
  };

  EXPECT_EQ(mapEnumString.size(), autopas::TuningStrategyOption::getOptionNames().size());

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

TEST(OptionTest, parseAcquisitionFunctionOptionsTest) {
  std::map<autopas::AcquisitionFunctionOption, std::string> mapEnumString = {
      {autopas::AcquisitionFunctionOption::upperConfidenceBound, "upperconfbound"},
      {autopas::AcquisitionFunctionOption::expectedImprovement, "expimprv"},
      {autopas::AcquisitionFunctionOption::mean, "mean"},
      {autopas::AcquisitionFunctionOption::probabilityOfImprovement, "probofimprv"},
      {autopas::AcquisitionFunctionOption::variance, "varianz"},
  };

  EXPECT_EQ(mapEnumString.size(), autopas::AcquisitionFunctionOption::getOptionNames().size());

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

// Generated tests for all option types
// parseOptionExact tests

TYPED_TEST_SUITE_P(OptionTest);

/**
 * tests parseOptionExact against the strings from getAllOptionNames.
 */
TYPED_TEST_P(OptionTest, parseExactOptionsTest) {
  auto mapOptionEnumString = TypeParam::getOptionNames();
  ASSERT_THAT(mapOptionEnumString, ::testing::SizeIs(::testing::Ge(1)));

  for (auto &[optionEnum, optionString] : mapOptionEnumString) {
    auto parsedOption = TypeParam::parseOptionExact(optionString);
    EXPECT_EQ(parsedOption, optionEnum);
  }
}

// to_string tests
/**
 * Test to to_string function for a given list of options
 */
TYPED_TEST_P(OptionTest, to_stringTest) {
  // good options
  for (auto &[optionEnum, optionString] : TypeParam::getOptionNames()) {
    EXPECT_EQ(optionEnum.to_string(), optionString);
  }

  // bad Options
  for (auto &op : {TypeParam()}) {
    EXPECT_THAT(op.to_string(), ::testing::HasSubstr("Unknown"));
  }
}

REGISTER_TYPED_TEST_SUITE_P(OptionTest, parseExactOptionsTest, to_stringTest);

// instantiate tests for all option types
using OptionTypes = ::testing::Types<autopas::AcquisitionFunctionOption, autopas::ContainerOption,
                                     autopas::DataLayoutOption, autopas::Newton3Option, autopas::SelectorStrategyOption,
                                     autopas::TraversalOption, autopas::TuningStrategyOption>;
INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, OptionTest, OptionTypes);
