/**
 * @file TraversalOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>
#include "autopas/options/Option.h"

namespace autopas {

class TraversalOption : public Option<TraversalOption> {
 public:

  /**
   * Possible choices for the cell pair traversal.
   */
  enum Value {
    c08 = 0,
    sliced = 1,
    c18 = 2,
    c01 = 3,
    directSumTraversal = 4,
    slicedVerlet = 5,
    c18Verlet = 6,
    c01Verlet = 7,
    c01Cuda = 8,
    verletTraversal = 9,
    c01CombinedSoA = 10,
    verletClusters = 11,
    c04 = 12,
    varVerletTraversalAsBuild = 13,
    verletClustersColoring = 14,
    c04SoA = 15,
  };

  TraversalOption() =  default;
  constexpr TraversalOption(Value option) : _value(option) {}
  constexpr operator Value() const { return _value; }
  explicit operator bool() = delete;

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   */
  static std::map<TraversalOption, std::string> getOptionNames() {
    return {
        {TraversalOption::c08, "c08"},
        {TraversalOption::sliced, "sliced"},
        {TraversalOption::c18, "c18"},
        {TraversalOption::c01, "c01"},
        {TraversalOption::directSumTraversal, "directSum"},
        {TraversalOption::slicedVerlet, "slicedVerlet"},
        {TraversalOption::c18Verlet, "verlet-c18"},
        {TraversalOption::c01Verlet, "verlet-c01"},
        {TraversalOption::c01Cuda, "c01Cuda"},
        {TraversalOption::verletTraversal, "verlet-lists"},
        {TraversalOption::c01CombinedSoA, "c01-combined-SoA"},
        {TraversalOption::verletClusters, "verlet-clusters"},
        {TraversalOption::c04, "c04"},
        {TraversalOption::varVerletTraversalAsBuild, "var-verlet-lists-as-build"},
        {TraversalOption::verletClustersColoring, "verlet-clusters-coloring"},
        {TraversalOption::c04SoA, "c04SoA"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
