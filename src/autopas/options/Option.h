/**
 * @file Option.h
 * @author F. Gratl
 * @date 10/14/19
 */

#pragma once

#include <map>
#include <set>

namespace autopas {

template <typename actualOption>
class Option {
 public:

  constexpr Option() = default;

  /**
   * Provides a way to iterate over the possible options.
   */
  static std::set<actualOption> getAllOptions() {
    std::set<actualOption> retSet;
    std::for_each(actualOption::getOptionNames().begin(), actualOption::getOptionNames().end(),
                   [&retSet](auto pairOpStr) { retSet.insert(pairOpStr.first); });
    return retSet;
  };

  std::string to_string() {
    auto &actualThis = *static_cast<actualOption *>(this);
    auto match = actualOption::getOptionsNames().find(actualThis);
    if (match == actualOption::getOptionsNames().end()) {
      return "Unknown TraversalOption (" + std::to_string(actualThis) + ")";
    } else {
      return match->second;
    }
  }
};
}  // namespace autopas
