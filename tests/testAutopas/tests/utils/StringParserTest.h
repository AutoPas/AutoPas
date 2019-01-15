/**
 * @file StringParserTest.h
 * @author F. Gratl
 * @date 1/15/19
 */

#pragma once

#include <gmock/gmock-matchers.h>
#include "AutoPasTestBase.h"
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"
#include "autopas/utils/StringParser.h"

class StringParserTest : public AutoPasTestBase {};

template <class T>
void testParseMultiple(std::vector<T> &allOptions, std::string &&optionsString,
                       std::function<std::vector<T>(std::string &)> &&parseFun) {
  std::sort(allOptions.begin(), allOptions.end());

  auto parsedOptions = parseFun(optionsString);

  std::sort(parsedOptions.begin(), parsedOptions.end());

  EXPECT_THAT(parsedOptions, ::testing::ContainerEq(allOptions));
}

template <class T>
void testParseSingle(std::vector<T> &allOptions, std::vector<std::string> &&optionsStrings,
                     std::function<T(std::string &)> parseFun) {
  std::sort(allOptions.begin(), allOptions.end());

  ASSERT_EQ(allOptions.size(), optionsStrings.size()) << "Not all options tested!";

  std::vector<T> parsedOptions;

  for (auto &string : optionsStrings) {
    parsedOptions.push_back(parseFun(string));
  }

  ASSERT_THAT(parsedOptions, ::testing::UnorderedElementsAreArray(allOptions));
}