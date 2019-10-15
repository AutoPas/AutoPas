/**
 * @file DataLayoutOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the traversal choices.
 */
class DataLayoutOption : public Option<DataLayoutOption> {
 public:
  /**
   * Possible choices for the particle data layout.
   */
  enum Value { aos, soa, cuda };

  DataLayoutOption() = default;
  constexpr DataLayoutOption(Value option) : _value(option) {}
  constexpr operator Value() const { return _value; }
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   */
  static std::map<DataLayoutOption, std::string> getOptionNames() {
    return {
        {DataLayoutOption::aos, "Array-of-Structures"},
        {DataLayoutOption::soa, "Structure-of-Arrays"},
        {DataLayoutOption::cuda, "Structure-of-Arrays on CUDA capable device"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
