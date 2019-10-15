/**
 * @file Option.h
 * @author F. Gratl
 * @date 10/14/19
 */

#pragma once

#include <map>
#include <set>
#include "autopas/utils/StringUtils.h"

namespace autopas {

template <typename actualOption>
class Option {
 public:
  constexpr Option() = default;

  /**
   * Provides a way to iterate over the possible options.
   * @return set of all possible values of this option type.
   */
  static std::set<actualOption> getAllOptions() {
    std::set<actualOption> retSet;
    std::for_each(actualOption::getOptionNames().begin(), actualOption::getOptionNames().end(),
                  [&retSet](auto pairOpStr) { retSet.insert(pairOpStr.first); });
    return retSet;
  };

  /**
   * Converts an Option object to its respective string representation.
   * @return The string representation or "Unknown Option (<IntValue>)".
   */
  std::string to_string() const {
    auto &actualThis = *static_cast<const actualOption *>(this);
    auto match = actualOption::getOptionNames().find(actualThis);
    if (match == actualOption::getOptionNames().end()) {
      return "Unknown Option (" + std::to_string(actualThis) + ")";
    } else {
      return match->second;
    }
  }

  /**
   * Converts a string of options to a set of enums. For best results, the options are expected to be lower case.
   *
   * Allowed delimiters can be found in autopas::utils::StringUtils::delimiters.
   * Possible options can be found in getAllOptions().
   *
   * This function uses the Needleman-Wunsch algorithm to find the closest matching options.
   * If an option is ambiguous an execption is thrown.
   *
   * @param optionsString String containing traversal options.
   * @return Set of option enums. If no valid option was found the empty set is returned.
   */
  static std::set<actualOption> parseOptions(const std::string &optionsString) {
    std::set<actualOption> optionsSet;

    auto needles = autopas::utils::StringUtils::tokenize(optionsString, autopas::utils::StringUtils::delimiters);
    std::vector<std::string> haystack;

    std::transform(actualOption::getAllOptions().begin(), actualOption::getAllOptions().end(),
                   std::back_inserter(haystack), [](auto traversalOption) { return traversalOption.to_string(); });

    // assure all string representations are lower case
    std::map<actualOption, std::string> allOptionNames;
    for (auto &pairEnumString : actualOption::getOptionNames()) {
      auto s = pairEnumString.second;
      std::transform(s.begin(), s.end(), s.begin(), ::tolower);
      allOptionNames.emplace(pairEnumString.first, s);
    }

    for (auto &needle : needles) {
      auto matchingString = autopas::utils::StringUtils::matchStrings(haystack, needle);
      for (auto &pairEnumString : allOptionNames) {
        if (pairEnumString.second == matchingString) {
          optionsSet.insert(pairEnumString.first);
        }
      }
    }

    return optionsSet;
  }
};
}  // namespace autopas
