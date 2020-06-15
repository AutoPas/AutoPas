/**
 * @file ContainerOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
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
    verletClusterCells = 6,
    referenceLinkedCells = 7,
  };

  /**
   * Constructor.
   */
  ContainerOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr ContainerOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * No cast to bool.
   * @return
   */
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of ContainerOption.
   * @return map option -> string representation
   */
  static std::map<ContainerOption, std::string> getOptionNames() {
    return {
        {ContainerOption::directSum, "DirectSum"},
        {ContainerOption::linkedCells, "LinkedCells"},
        {ContainerOption::referenceLinkedCells, "ReferenceLinkedCells"},
        {ContainerOption::verletLists, "VerletLists"},
        {ContainerOption::verletListsCells, "VerletListsCells"},
        {ContainerOption::verletClusterLists, "VerletClusterLists"},
        {ContainerOption::varVerletListsAsBuild, "VarVerletListsAsBuild"},
        {ContainerOption::verletClusterCells, "VerletClusterCells"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
