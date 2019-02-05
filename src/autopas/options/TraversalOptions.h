/**
 * @file TraversalOptions.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <vector>

namespace autopas {

/**
 * Possible choices for the cell pair traversal.
 */
enum TraversalOptions {
  c08 = 0,
  sliced = 1,
  c18 = 2,
  c01 = 3,
  directSumTraversal = 4,
  slicedVerlet = 5,
  c18Verlet = 6,
  c01Verlet = 7,
  dummyTraversal = 666,
  directSumTraversalCuda = 8,
};

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static const std::vector<TraversalOptions> allTraversalOptions = {
    TraversalOptions::c08,
    TraversalOptions::sliced,
    TraversalOptions::c18,
    TraversalOptions::c01,
    TraversalOptions::directSumTraversal,
    TraversalOptions::slicedVerlet,
    TraversalOptions::c18Verlet,
    TraversalOptions::c01Verlet,
};

}  // namespace autopas
