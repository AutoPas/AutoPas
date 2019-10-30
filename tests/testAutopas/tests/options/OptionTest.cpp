/**
 * @file OptionTest.cpp
 * @author F. Gratl
 * @date 10/29/19
 */

#include "OptionTest.h"
#include <autopas/options/AcquisitionFunctionOption.h>
#include <autopas/options/Newton3Option.h>
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "tests/utils/StringUtilsTest.h"

// parseOptions tests
// these tests shall not (yet :) be generated since we here want to pass strings that do not match exactly.

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

TEST(OptionTest, parseSelectorOptionsTest) {
  testParseSingle<autopas::SelectorStrategyOption>(autopas::SelectorStrategyOption::getAllOptions(),
                                                   {"absolute", "median", "mean"},
                                                   autopas::SelectorStrategyOption::parseOptions);
}

TEST(OptionTest, parseTuningStrategyOptionsTest) {
  testParseSingle<autopas::TuningStrategyOption>(autopas::TuningStrategyOption::getAllOptions(),
                                                 {"full-search", "bayesian-search"},
                                                 autopas::TuningStrategyOption::parseOptions);
}

// Generated tests for all option types
// parseOptionExact tests

TYPED_TEST_SUITE_P(OptionTest);

TYPED_TEST_P(OptionTest, parseExactOptionsTest) {
  std::vector<std::string> allOptionNames;
  allOptionNames.reserve(TypeParam::getAllOptions().size());

  for (auto [_, optionName] : TypeParam::getOptionNames()) {
    allOptionNames.push_back(optionName);
  }
  testParseExact<TypeParam>(TypeParam::getAllOptions(), allOptionNames, TypeParam::parseOptionExact);
}

// to_string tests
TYPED_TEST_P(OptionTest, to_stringTest) {
  testToString(TypeParam::getAllOptions(), {TypeParam()});
}

REGISTER_TYPED_TEST_SUITE_P(OptionTest, parseExactOptionsTest, to_stringTest);

// instantiate tests for all option types
typedef ::testing::Types<autopas::AcquisitionFunctionOption, autopas::ContainerOption, autopas::DataLayoutOption,
                         autopas::Newton3Option, autopas::SelectorStrategyOption, autopas::TraversalOption,
                         autopas::TuningStrategyOption>
    OptionTypes;
INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, OptionTest, OptionTypes);
