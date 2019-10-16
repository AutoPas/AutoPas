/**
 * @file StringUtilsTest.cpp
 * @author F. Gratl
 * @date 1/15/19
 */

#include "StringUtilsTest.h"
#include <autopas/options/ContainerOption.h>
#include <autopas/options/DataLayoutOption.h>
#include <autopas/options/SelectorStrategyOption.h>
#include <autopas/options/TraversalOption.h>
#include <autopas/options/TuningStrategyOption.h>

TEST(StringUtilsTest, parseTraversalOptionsTest) {
  testParseMultiple<autopas::TraversalOption>(
      autopas::TraversalOption::getAllOptions(),
      "c01, c04, c08, c18, c04s, direct, sliced v01, c18verlet, verlet-sliced, "
      "cuda-c01, verlet-lists, c01-combined, verlet-clusters, var-verlet-lists-as-build, verlet-clusters-coloring, "
      "verlet-cluster-cells",
      autopas::TraversalOption::parseOptions);
}

TEST(StringUtilsTest, parseContainerOptionsTest) {
  testParseMultiple<autopas::ContainerOption>(
      autopas::ContainerOption::getAllOptions(),
      "directSum, linkedCells, verletLists, verletLists-cells, vclusterlists, varVerletListsAsBuild, vclustercells",
      autopas::ContainerOption::parseOptions);
}

TEST(StringUtilsTest, parseDataLayoutOptionsTest) {
#if defined(AUTOPAS_CUDA)
  auto options = "cuda, soa, aos";
#else
  auto options = "soa, aos";
#endif
  testParseMultiple<autopas::DataLayoutOption>(autopas::DataLayoutOption::getAllOptions(), options,
                                               autopas::DataLayoutOption::parseOptions);
}

TEST(StringUtilsTest, parseDoublesTest) {
  testParseMultiple<double>({1., 1.5, 2., 3., 20.}, "1.,1.5, 2,3.00,2e1", autopas::utils::StringUtils::parseDoubles);
}

TEST(StringUtilsTest, parseNumberSetTest) {
  EXPECT_EQ(autopas::utils::StringUtils::parseNumberSet("1.,1.5, 2,3.00,2e1")->getAll(),
            std::set<double>({1., 1.5, 2., 3., 20.}));

  auto numberSet = autopas::utils::StringUtils::parseNumberSet("[1.,2e1]");
  auto *numberInterval = dynamic_cast<autopas::NumberInterval<double> *>(numberSet.get());
  EXPECT_NE(numberInterval, nullptr);
  if (numberInterval) {
    EXPECT_EQ(numberInterval->getMin(), 1.);
    EXPECT_EQ(numberInterval->getMax(), 2e1);
  }
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
