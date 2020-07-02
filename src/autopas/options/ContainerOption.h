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
    directSum,
    linkedCells,
    verletLists,
    verletListsCells,
    verletClusterLists,
    varVerletListsAsBuild,
    verletClusterCells,
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
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<ContainerOption> getDiscouragedOptions() { return {Value::directSum, Value::verletClusterCells}; }

  /**
   * Provides a way to iterate over the possible choices of ContainerOption.
   * @return map option -> string representation
   */
  static std::map<ContainerOption, std::string> getOptionNames() {
    return {
        {ContainerOption::directSum, "DirectSum"},
        {ContainerOption::linkedCells, "LinkedCells"},
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
