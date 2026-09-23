/**
 * @file YamlParserTest.cpp
 * @author S. Newcome
 * @date 26.06.2026
 */

#include "YamlParserTest.h"

#include <yaml-cpp/yaml.h>

#include "src/configuration/YamlParser.h"
#include "src/domainDecomposition/LoadBalancerOption.h"

using MDFlexParser::YamlParser::parseSequenceOneElementExpected;

/**
 * A scalar node holding a single value is returned unchanged.
 */
TEST_F(YamlParserTest, scalarReturnsValue) {
  const auto node = YAML::Load("InvertedPressure");
  EXPECT_EQ(parseSequenceOneElementExpected(node, "err"), "InvertedPressure");
}

/**
 * A sequence holding exactly one value returns that value (without sequence decoration).
 */
TEST_F(YamlParserTest, singleElementSequenceReturnsValue) {
  const auto node = YAML::Load("[InvertedPressure]");
  EXPECT_EQ(parseSequenceOneElementExpected(node, "err"), "InvertedPressure");
}

/**
 * A sequence with more than one element violates the single-element expectation and throws.
 */
TEST_F(YamlParserTest, multiElementSequenceThrows) {
  const auto node = YAML::Load("[InvertedPressure, None]");
  EXPECT_THROW(parseSequenceOneElementExpected(node, "err"), std::runtime_error);
}

/**
 * By default "all" is rejected: parsed as an option it would expand into every option, violating the
 * single-element expectation.
 */
TEST_F(YamlParserTest, allThrowsByDefault) {
  EXPECT_THROW(parseSequenceOneElementExpected(YAML::Load("all"), "err"), std::runtime_error);
  EXPECT_THROW(parseSequenceOneElementExpected(YAML::Load("[all]"), "err"), std::runtime_error);
}

/**
 * With allThrowsError disabled, "all" is returned verbatim instead of throwing. This is the path the load
 * balancer relies on, where "all" is meant as the single option "A Load-balancing Library (ALL)".
 */
TEST_F(YamlParserTest, allReturnedWhenNotThrowing) {
  EXPECT_EQ(parseSequenceOneElementExpected(YAML::Load("all"), "err", false), "all");
  EXPECT_EQ(parseSequenceOneElementExpected(YAML::Load("[all]"), "err", false), "all");
}

/**
 * End-to-end check of the load balancer special case: the user's "all" is upper-cased to "ALL" and parsed into
 * LoadBalancerOption::all rather than being rejected or expanded into every option.
 */
TEST_F(YamlParserTest, loadBalancerAllParsedAsALL) {
  auto loadBalancerString = parseSequenceOneElementExpected(YAML::Load("[all]"), "err", false);
  if (loadBalancerString == "all") {
    loadBalancerString = "ALL";
  }
  const auto parsedOptions = LoadBalancerOption::parseOptions(loadBalancerString);
  ASSERT_EQ(parsedOptions.size(), 1);
  EXPECT_EQ(*parsedOptions.begin(), LoadBalancerOption::all);
}