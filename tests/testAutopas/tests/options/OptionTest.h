/**
 * @file OptionTest.h
 * @author F. Gratl
 * @date 29.10.2019
 */

#pragma once

#include <gmock/gmock-matchers.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ArrayUtils.h"

namespace OptionTest {

/**
 * Tests for all Options derived from autopas::Option.
 * @tparam T Template needed for Type-Parameterized Tests
 */
template <typename T>
class OptionTest : public AutoPasTestBase {};

/**
 * For each map entry it is tested if the string can be parsed to the respective enum with parseOptions.
 * @tparam T Type of the object resulting form parsing.
 * @param mapOptionString mapping of enums to strings.
 */
template <class T>
void testParseOptionsIndividually(const std::map<T, std::string> &mapOptionString) {
  for (auto &[optionEnum, optionString] : mapOptionString) {
    auto parsedOptions = T::parseOptions(optionString);
    EXPECT_THAT(parsedOptions, ::testing::SizeIs(1))
        << "Option " << optionEnum.to_string()
        << " was parsed ambiguously: " << autopas::utils::ArrayUtils::to_string(parsedOptions);
    EXPECT_EQ(*(parsedOptions.begin()), optionEnum)
        << "Option " << optionEnum.to_string() << " was not correctly parsed!";
  }
}

/**
 * Tests if parseOptions can parse all options when all strings of the map are combined to one.
 * @tparam T Type of the object resulting form parsing.
 * @param mapOptionString mapping of enums to strings.
 */
template <class T>
void testParseOptionsCombined(const std::map<T, std::string> &mapOptionString) {
  std::ostringstream allOptionsStringStream;

  // merge all strings in one separated by a legal delimiter. We don't care for the trailing comma
  for (auto &[_, optionString] : mapOptionString) {
    allOptionsStringStream << optionString << ", ";
  }

  auto parsedOptions = T::parseOptions(allOptionsStringStream.str());
  EXPECT_EQ(parsedOptions.size(), mapOptionString.size())
      << "Incorrect number of options parsed! Following options were found: "
      << autopas::utils::ArrayUtils::to_string(parsedOptions);
}

}  // end namespace OptionTest
