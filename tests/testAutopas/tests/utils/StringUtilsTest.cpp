/**
 * @file StringUtilsTest.cpp
 * @author F. Gratl
 * @date 1/15/19
 */

#include "StringUtilsTest.h"

TEST(StringUtilsTest, parseTraversalOptionsTest) {
  testParseMultiple<autopas::TraversalOption>(
      autopas::allTraversalOptions,
      "c01, c08, c18, direct; sliced v01, c18verlet, verlet-sliced, cuda-c01, verlet-lists",
      autopas::utils::StringUtils::parseTraversalOptions);
}

TEST(StringUtilsTest, parseContainerOptionsTest) {
  testParseMultiple<autopas::ContainerOption>(autopas::allContainerOptions,
                                              "directSum, linkedCells, verletLists, verlet-cells, vcluster",
                                              autopas::utils::StringUtils::parseContainerOptions);
}

TEST(StringUtilsTest, parseDataLayoutOptionsTest) {
#if defined(AUTOPAS_CUDA)
  testParseMultiple<autopas::DataLayoutOption>(autopas::allDataLayoutOptions, "cuda, soa, aos",
                                               autopas::utils::StringUtils::parseDataLayout);
#else
  testParseMultiple<autopas::DataLayoutOption>(autopas::allDataLayoutOptions, "soa, aos",
                                               autopas::utils::StringUtils::parseDataLayout);
#endif
}

TEST(StringUtilsTest, parseDoublesTest) {
  testParseMultiple<double>({1., 1.5, 2., 3., 20.}, "1.,1.5, 2,3.00,2e1", autopas::utils::StringUtils::parseDoubles);
}

TEST(StringUtilsTest, parseNumberSetTest) {
  EXPECT_EQ(autopas::utils::StringUtils::parseNumberSet("1.,1.5, 2,3.00,2e1")->getAll(),
            std::set<double>({1., 1.5, 2., 3., 20.}));

  auto numberSet = autopas::utils::StringUtils::parseNumberSet("[1.,2e1]");
  auto* numberInterval = dynamic_cast<autopas::NumberInterval<double>*>(numberSet.get());
  EXPECT_NE(numberInterval, nullptr);
  if (numberInterval) {
    EXPECT_EQ(numberInterval->getMin(), 1.);
    EXPECT_EQ(numberInterval->getMax(), 2e1);
  }
}

TEST(StringUtilsTest, parseSelectorOptionsTest) {
  testParseSingle<autopas::SelectorStrategyOption>(autopas::allSelectorStrategies, {"absolute", "median", "mean"},
                                                   autopas::utils::StringUtils::parseSelectorStrategy);
}

TEST(StringUtilsTest, parseTuningStrategyOptionsTest) {
  testParseSingle<autopas::TuningStrategyOption>(autopas::allTuningStrategyOptions, {"full-search"},
                                                 autopas::utils::StringUtils::parseTuningStrategyOption);
}

TEST(StringUtilsTest, to_stringDataLayoutTest) {
  testToString(autopas::allDataLayoutOptions, {autopas::DataLayoutOption(-1)});
}

TEST(StringUtilsTest, to_stringSelectorStrategiesTest) {
  testToString(autopas::allSelectorStrategies, {autopas::SelectorStrategyOption(-1)});
}

TEST(StringUtilsTest, to_stringContainerOptionsTest) {
  testToString(autopas::allContainerOptions, {autopas::ContainerOption(-1)});
}

TEST(StringUtilsTest, to_stringTraversalOptionsTest) {
  testToString(autopas::allTraversalOptions, {autopas::TraversalOption(-1)});
}

TEST(StringUtilsTest, to_stringTuningStrategyOptionsTest) {
  testToString(autopas::allTuningStrategyOptions, {autopas::TuningStrategyOption(-1)});
}
