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
#include "autopas/options/OpenMPKindOption.h"
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
  // collect all option names
  std::map<autopas::TraversalOption, std::string> mapEnumString = autopas::TraversalOption::getOptionNames();

  // alter all strings
  std::transform(mapEnumString.begin(), mapEnumString.end(), std::inserter(mapEnumString, mapEnumString.end()),
                 [](auto &pair) {
                   auto &[option, str] = pair;
                   // to lower
                   std::transform(str.begin(), str.end(), str.begin(), ::tolower);
                   // remove all underscores
                   str.erase(std::remove(str.begin(), str.end(), '_'), str.end());
                   return std::make_pair(option, str);
                 });

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

TEST(OptionTest, parseContainerOptionsTest) {
  std::map<autopas::ContainerOption, std::string> mapEnumString = {
      {autopas::ContainerOption::directSum, "directSum"},
      {autopas::ContainerOption::linkedCells, "linkedCells"},
      {autopas::ContainerOption::varVerletListsAsBuild, "varVerletListsAsBuild"},
      {autopas::ContainerOption::verletClusterLists, "vclusterlists"},
      {autopas::ContainerOption::verletLists, "verletLists"},
      {autopas::ContainerOption::verletListsCells, "verletLists-cells"},
      {autopas::ContainerOption::linkedCellsReferences, "linkedCellsreferenc"},
      {autopas::ContainerOption::pairwiseVerletLists, "pairwiseVerlet"},
      {autopas::ContainerOption::octree, "octree"}};

  EXPECT_EQ(mapEnumString.size(), autopas::ContainerOption::getOptionNames().size());

  testParseOptionsIndividually(mapEnumString);
  testParseOptionsCombined(mapEnumString);
}

TEST(OptionTest, parseDataLayoutOptionsTest) {
  std::map<autopas::DataLayoutOption, std::string> mapEnumString = {
      {autopas::DataLayoutOption::aos, "aos"},
      {autopas::DataLayoutOption::soa, "soa"},
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
      {autopas::TuningStrategyOption::activeHarmony, "harmony"},
      {autopas::TuningStrategyOption::bayesianSearch, "bayesian"},
      {autopas::TuningStrategyOption::bayesianClusterSearch, "bayesian-cluster"},
      {autopas::TuningStrategyOption::fullSearch, "full"},
      {autopas::TuningStrategyOption::mpiDivideAndConquer, "divide&conquer"},
      {autopas::TuningStrategyOption::predictiveTuning, "predictive"},
      {autopas::TuningStrategyOption::randomSearch, "random"},
      {autopas::TuningStrategyOption::ruleBasedTuning, "rule-based"},
      {autopas::TuningStrategyOption::fuzzyTuning, "fuzzy"},
      {autopas::TuningStrategyOption::slowConfigFilter, "slow-filter"},
      {autopas::TuningStrategyOption::sortByName, "sortbyname"},
      {autopas::TuningStrategyOption::tuningStrategyLogger, "tuningstratLogger"},
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

TEST(OptionTest, parseOpenMPKindOptionsTest) {
  std::map<autopas::OpenMPKindOption, std::string> mapEnumString = {
      // Standard OpenMP's scheduling kinds:
      {autopas::OpenMPKindOption::omp_auto, "auto"},
      {autopas::OpenMPKindOption::omp_dynamic, "dynamic"},
      {autopas::OpenMPKindOption::omp_guided, "guided"},
      {autopas::OpenMPKindOption::omp_runtime, "runtime"},
      {autopas::OpenMPKindOption::omp_static, "static"},

      // Auto4OMP's automated selection methods:
      {autopas::OpenMPKindOption::auto4omp_randomsel, "randomSel"},
      {autopas::OpenMPKindOption::auto4omp_exhaustivesel, "exhaustiveSel"},
      {autopas::OpenMPKindOption::auto4omp_binarySearch, "binarySearch"},
      {autopas::OpenMPKindOption::auto4omp_expertsel, "expertSel"},

#ifdef AUTOPAS_USE_LB4OMP
      // LB4OMP's scheduling techniques (beware, technique names in LB4OMP's README are outdated):
      {autopas::OpenMPKindOption::lb4omp_profiling, "profiling"},  // Profiling
      {autopas::OpenMPKindOption::lb4omp_fsc, "fsc"},              // Fixed Size Chunk
      {autopas::OpenMPKindOption::lb4omp_mfsc, "mfsc"},            // Modified Fixed Size Chunk
      {autopas::OpenMPKindOption::lb4omp_tap, "tap"},              // Tapering
      {autopas::OpenMPKindOption::lb4omp_fac, "fac"},              // Factoring
      {autopas::OpenMPKindOption::lb4omp_faca, "faca"},            // Improved Factoring
      {autopas::OpenMPKindOption::lb4omp_bold, "bold"},            // Bold
      {autopas::OpenMPKindOption::lb4omp_fac2, "fac2"},            // Practical Factoring
      {autopas::OpenMPKindOption::lb4omp_wf, "wf"},                // Weighted Factoring
      {autopas::OpenMPKindOption::lb4omp_af, "af"},                // Adaptive Factoring
      {autopas::OpenMPKindOption::lb4omp_awf, "awf"},              // Adaptive Weighted Factoring
      {autopas::OpenMPKindOption::lb4omp_tfss, "tfss"},            // Trapezoid Factoring Self Scheduling
      {autopas::OpenMPKindOption::lb4omp_fiss, "fiss"},            // Fixed Increase Self Scheduling
      {autopas::OpenMPKindOption::lb4omp_viss, "viss"},            // Variable Increase Self Scheduling
      {autopas::OpenMPKindOption::lb4omp_rnd, "rnd"},              // Random

      // LB4OMP's scheduling techniques used by Auto4OMP (in addition to the standard scheduling kinds):
      {autopas::OpenMPKindOption::lb4omp_trapezoidal, "trapezoidal"},    // Trapezoid Self Scheduling
      {autopas::OpenMPKindOption::lb4omp_fac2a, "fac2a"},                // Improved Practical Factoring
      {autopas::OpenMPKindOption::lb4omp_static_steal, "static_steal"},  // Static with Steal enabled
      {autopas::OpenMPKindOption::lb4omp_awf_b, "awf_b"},                // Adaptive Weighted Factoring Variant B
      {autopas::OpenMPKindOption::lb4omp_awf_c, "awf_c"},                // Adaptive Weighted Factoring Variant C
      {autopas::OpenMPKindOption::lb4omp_awf_d, "awf_d"},                // Adaptive Weighted Factoring Variant D
      {autopas::OpenMPKindOption::lb4omp_awf_e, "awf_e"},                // Adaptive Weighted Factoring Variant E
      {autopas::OpenMPKindOption::lb4omp_af_a, "af_a"},                  // Improved Adaptive Factoring
#endif
  };

  EXPECT_EQ(mapEnumString.size(), autopas::OpenMPKindOption::getOptionNames().size());

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
using OptionTypes =
    ::testing::Types<autopas::AcquisitionFunctionOption, autopas::ContainerOption, autopas::DataLayoutOption,
                     autopas::Newton3Option, autopas::SelectorStrategyOption, autopas::TraversalOption,
                     autopas::TuningStrategyOption, autopas::OpenMPKindOption>;
INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, OptionTest, OptionTypes);
