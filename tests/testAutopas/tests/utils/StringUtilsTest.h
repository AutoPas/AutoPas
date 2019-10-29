/**
 * @file StringUtilsTest.h
 * @author F. Gratl
 * @date 1/15/19
 */

#pragma once

#include <gmock/gmock-matchers.h>
#include <functional>
#include "AutoPasTestBase.h"
#include "autopas/utils/StringUtils.h"

class StringUtilsTest : public AutoPasTestBase {};

/**
 * Tests a parsing function which takes a string and returns a set of values.
 * @tparam T Type of the object resulting form parsing.
 * @param allOptions Set of all expected options.
 * @param optionsString String that shall be parsed.
 * @param parseFun Parsing Function.
 */
template <class T>
void testParseMultiple(const std::set<T> &allOptions, const std::string &optionsString,
                       std::function<std::set<T>(const std::string &)> &&parseFun) {
  auto parsedOptions = parseFun(optionsString);

  EXPECT_THAT(parsedOptions, ::testing::ElementsAreArray(allOptions));
}
