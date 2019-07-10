/**
 * @file TraversalOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

namespace autopas {

/**
 * Possible choices for the cell pair traversal.
 */
enum TraversalOption {
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
};

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static const std::set<TraversalOption> allTraversalOptions = {
    TraversalOption::c08,
    TraversalOption::sliced,
    TraversalOption::c18,
    TraversalOption::c01,
    TraversalOption::directSumTraversal,
    TraversalOption::slicedVerlet,
    TraversalOption::c18Verlet,
    TraversalOption::c01Verlet,
    TraversalOption::c01Cuda,
    TraversalOption::verletTraversal,
    TraversalOption::c01CombinedSoA,
    TraversalOption::verletClusters,
    TraversalOption::c04,
};

}  // namespace autopas
