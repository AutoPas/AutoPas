/**
 * @file OptionTest.cpp
 * @author F. Gratl
 * @date 10/29/19
 */

#include "OptionTest.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "tests/utils/StringUtilsTest.h"

// parseOptions tests

TEST(OptionTest, parseTraversalOptionsTest) {
  testParseMultiple<autopas::TraversalOption>(
      autopas::TraversalOption::getAllOptions(),
      "c01, c04, c08, c18, c04-soa, direct, slicedv01, verletc18, verlec01, verlet-sliced, "
      "cudac01, verletlists, c01-combined, verlet-clusters, var-verlet-lists-as-build, verlet-clusters-coloring, "
      "verlet-cluster-cells",
      autopas::TraversalOption::parseOptions);
}

TEST(OptionTest, parseContainerOptionsTest) {
  testParseMultiple<autopas::ContainerOption>(
      autopas::ContainerOption::getAllOptions(),
      "directSum, linkedCells, verletLists, verletLists-cells, vclusterlists, varVerletListsAsBuild, vclustercells",
      autopas::ContainerOption::parseOptions);
}

TEST(OptionTest, parseDataLayoutOptionsTest) {
#if defined(AUTOPAS_CUDA)
  auto options = "cuda, soa, aos";
#else
  auto options = "soa, aos";
#endif
  testParseMultiple<autopas::DataLayoutOption>(autopas::DataLayoutOption::getAllOptions(), options,
                                               autopas::DataLayoutOption::parseOptions);
}

TEST(StringUtilsTest, parseSelectorOptionsTest) {
  testParseSingle<autopas::SelectorStrategyOption>(autopas::SelectorStrategyOption::getAllOptions(),
                                                   {"absolute", "median", "mean"},
                                                   autopas::SelectorStrategyOption::parseOptions);
}

TEST(StringUtilsTest, parseTuningStrategyOptionsTest) {
  testParseSingle<autopas::TuningStrategyOption>(autopas::TuningStrategyOption::getAllOptions(),
                                                 {"full-search", "bayesian-search"},
                                                 autopas::TuningStrategyOption::parseOptions);
}

// parseOptionExact tests

TEST(OptionTest, parseExactTraversalOptionsTest) {
  std::vector<std::string> allOptionNames;
  allOptionNames.reserve(autopas::TraversalOption::getAllOptions().size());

  for (auto [_, optionName] : autopas::TraversalOption::getOptionNames()) {
    allOptionNames.push_back(optionName);
  }
  testParseExact<autopas::TraversalOption>(autopas::TraversalOption::getAllOptions(), allOptionNames,
                                           autopas::TraversalOption::parseOptionExact<false>);
}

// TEST(OptionTest, parseContainerOptionsTest) {
//  testParseMultiple<autopas::ContainerOption>(
//      autopas::ContainerOption::getAllOptions(),
//      "directSum, linkedCells, verletLists, verletLists-cells, vclusterlists, varVerletListsAsBuild, vclustercells",
//      autopas::ContainerOption::parseOptions);
//}
//
// TEST(OptionTest, parseDataLayoutOptionsTest) {
//#if defined(AUTOPAS_CUDA)
//  auto options = "cuda, soa, aos";
//#else
//  auto options = "soa, aos";
//#endif
//  testParseMultiple<autopas::DataLayoutOption>(autopas::DataLayoutOption::getAllOptions(), options,
//                                               autopas::DataLayoutOption::parseOptions);
//}
//
// TEST(StringUtilsTest, parseSelectorOptionsTest) {
//  testParseSingle<autopas::SelectorStrategyOption>(autopas::SelectorStrategyOption::getAllOptions(),
//                                                   {"absolute", "median", "mean"},
//                                                   autopas::SelectorStrategyOption::parseOptions);
//}
//
// TEST(StringUtilsTest, parseTuningStrategyOptionsTest) {
//  testParseSingle<autopas::TuningStrategyOption>(autopas::TuningStrategyOption::getAllOptions(),
//                                                 {"full-search", "bayesian-search"},
//                                                 autopas::TuningStrategyOption::parseOptions);
//}

// to_string tests

TEST(StringUtilsTest, to_stringDataLayoutTest) {
  testToString(autopas::DataLayoutOption::getAllOptions(), {autopas::DataLayoutOption()});
}

TEST(StringUtilsTest, to_stringSelectorStrategiesTest) {
  testToString(autopas::SelectorStrategyOption::getAllOptions(), {autopas::SelectorStrategyOption()});
}

TEST(StringUtilsTest, to_stringContainerOptionsTest) {
  testToString(autopas::ContainerOption::getAllOptions(), {autopas::ContainerOption()});
}

TEST(StringUtilsTest, to_stringTraversalOptionsTest) {
  // testing for bad options does not make sense anymore
  testToString(autopas::TraversalOption::getAllOptions(), {});
}

TEST(StringUtilsTest, to_stringTuningStrategyOptionsTest) {
  testToString(autopas::TuningStrategyOption::getAllOptions(), {autopas::TuningStrategyOption()});
}