/**
 * @file OptionTest.h
 * @author F. Gratl
 * @date 10/29/19
 */

#pragma once

#include <gmock/gmock-matchers.h>

#include "AutoPasTestBase.h"

/**
 * Tests for all Options derived from autopas::Option.
 * @tparam T Template needed for Type-Parameterized Tests
 */
template <typename T>
class OptionTest : public AutoPasTestBase {};

/**
 * Tests a parsing function which takes a string and returns a value in a set of size one.
 * @tparam T Type of the object resulting form parsing.
 * @param allOptions Set of all expected options.
 * @param optionsStrings Vector of strings that shall be parsed.
 * @param parseFun Parsing function.
 */
template <class T>
void testParseSingle(const std::set<T> &allOptions, const std::vector<std::string> &optionsStrings,
                     std::function<std::set<T>(const std::string &)> &&parseFun) {
  ASSERT_EQ(allOptions.size(), optionsStrings.size()) << "Not all options tested!";

  std::set<T> allParsedOptions;

  for (auto &string : optionsStrings) {
    auto parsedOptions = parseFun(string);
    EXPECT_THAT(parsedOptions, ::testing::SizeIs(1));
    allParsedOptions.insert(*parsedOptions.begin());
  }

  ASSERT_THAT(allParsedOptions, ::testing::ElementsAreArray(allOptions));
}
