/**
 * @file Option.h
 * @author F. Gratl
 * @date 10/14/19
 */

#pragma once

#include "autopas/utils/StringUtils.h"
#include <map>
#include <set>

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
  std::string to_string() {
    auto &actualThis = *static_cast<actualOption *>(this);
    auto match = actualOption::getOptionsNames().find(actualThis);
    if (match == actualOption::getOptionsNames().end()) {
      return "Unknown Option (" + std::to_string(actualThis) + ")";
    } else {
      return match->second;
    }
  }

  /**
   * Converts a string of options to a set of enums. The options are expected to be lower case.
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

    std::transform(actualOption::ge<tAllOptions().begin(), actualOption::getAllOptions().end(),
                   std::back_inserter(haystack),
                   [](auto traversalOption) { return actualOption::to_string(traversalOption); });

    for (auto &needle : needles) {
      auto matchingString = autopas::utils::StringUtils::matchStrings(haystack, needle);
      for (auto &pairEnumString : actualOption::getOptionNames()) {
        if (pairEnumString.second == matchingString) optionsSet.insert(pairEnumString.first);
      }
    }

    return optionsSet;
  }
};
}  // namespace autopas
