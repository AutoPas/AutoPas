/**
 * @file StringParserTest.h
 * @author F. Gratl
 * @date 1/15/19
 */

#pragma once

#include <gmock/gmock-matchers.h>
#include <functional>
#include "AutoPasTestBase.h"
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"
#include "autopas/utils/StringParser.h"

class StringParserTest : public AutoPasTestBase {};

/**
 * Tests a parsing function which takes a string and returns a vector of values.
 * @tparam T Type of the object resulting form parsing.
 * @param allOptions Vector of all expected options.
 * @param optionsString String that shall be parsed.
 * @param parseFun Parsing Function.
 */
template <class T>
void testParseMultiple(const std::vector<T> &allOptions, const std::string &optionsString,
                       std::function<std::vector<T>(const std::string &)> &&parseFun) {
  auto parsedOptions = parseFun(optionsString);

  std::sort(parsedOptions.begin(), parsedOptions.end());

  EXPECT_THAT(parsedOptions, ::testing::ContainerEq(allOptions));
}

/**
 * Tests a parsing function which takes a string and returns a value.
 * @tparam T Type of the object resulting form parsing.
 * @param allOptions Vector of all expected options.
 * @param optionsStrings Vector of strings that shall be parsed.
 * @param parseFun Parsing function.
 */
template <class T>
void testParseSingle(const std::vector<T> &allOptions, const std::vector<std::string> &optionsStrings,
                     std::function<T(const std::string &)> &&parseFun) {
  ASSERT_EQ(allOptions.size(), optionsStrings.size()) << "Not all options tested!";

  std::vector<T> parsedOptions;

  for (auto &string : optionsStrings) {
    parsedOptions.push_back(parseFun(string));
  }

  ASSERT_THAT(parsedOptions, ::testing::UnorderedElementsAreArray(allOptions));
}