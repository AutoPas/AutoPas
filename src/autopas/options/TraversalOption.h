/**
 * @file TraversalOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <vector>

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
  dummyTraversal = 666,
  directSumKokkosTraversal = 8,
};

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static const std::vector<TraversalOption> allTraversalOptions = {
    TraversalOption::c08,
    TraversalOption::sliced,
    TraversalOption::c18,
    TraversalOption::c01,
    TraversalOption::directSumTraversal,
    TraversalOption::slicedVerlet,
    TraversalOption::c18Verlet,
    TraversalOption::c01Verlet,
    TraversalOption::directSumKokkosTraversal,
};

}  // namespace autopas