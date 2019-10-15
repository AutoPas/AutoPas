/**
 * @file ContainerOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the container choices.
 */
class ContainerOption : public Option<ContainerOption> {
 public:
  /**
   * Possible choices for the particle container type.
   */
  enum Value {
    directSum = 0,
    linkedCells = 1,
    verletLists = 2,
    verletListsCells = 3,
    verletClusterLists = 4,
    varVerletListsAsBuild = 5,
  };

  ContainerOption() = default;
  constexpr ContainerOption(Value option) : _value(option) {}
  constexpr operator Value() const { return _value; }
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   */
  static std::map<ContainerOption, std::string> getOptionNames() {
    return {
        {ContainerOption::directSum, "DirectSum"},
        {ContainerOption::linkedCells, "LinkedCells"},
        {ContainerOption::verletLists, "VerletLists"},
        {ContainerOption::verletListsCells, "VerletListsCells"},
        {ContainerOption::verletClusterLists, "VerletClusterLists"},
        {ContainerOption::varVerletListsAsBuild, "VarVerletListsAsBuild"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas