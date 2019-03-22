/**
 * @file StringUtilsTest.cpp
 * @author F. Gratl
 * @date 1/15/19
 */

#include "StringUtilsTest.h"

TEST(StringUtilsTest, parseTraversalOptionsTest) {
  testParseMultiple<autopas::TraversalOption>(autopas::allTraversalOptions,
                                              "c01, c08, c18, direct; sliced v01, c18verlet, verlet-sliced, cuda-c01",
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

TEST(StringUtilsTest, parseSelectorOptionsTest) {
  testParseSingle<autopas::SelectorStrategy>(autopas::allSelectorStrategies, {"absolute", "median", "mean"},
                                             autopas::utils::StringUtils::parseSelectorStrategy);
}

TEST(StringUtilsTest, to_stringDataLayoutTest) {
  testToString(autopas::allDataLayoutOptions, {autopas::DataLayoutOption(-1)});
}

TEST(StringUtilsTest, to_stringSelectorStrategiesTest) {
  testToString(autopas::allSelectorStrategies, {autopas::SelectorStrategy(-1)});
}

TEST(StringUtilsTest, to_stringContainerOptionsTest) {
  testToString(autopas::allContainerOptions, {autopas::ContainerOption(-1)});
}

TEST(StringUtilsTest, to_stringTraversalOptionsTest) {
  testToString(autopas::allTraversalOptions, {autopas::TraversalOption(-1)});
}
