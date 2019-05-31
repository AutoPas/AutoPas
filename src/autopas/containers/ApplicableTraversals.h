/**
 * @file ApplicableTraversals.h
 * @author F. Gratl
 * @date 5/31/19
 */

#pragma once

#include <set>
#include "autopas/options/ContainerOption.h"
#include "autopas/options/TraversalOption.h"
namespace autopas {
namespace applicableTraversals {

/**
 * Lists all traversal options applicable for the Direct Sum container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allDSApplicableTraversals() {
  static const std::set<TraversalOption> v{TraversalOption::directSumTraversal};
  return v;
}

/**
 * Lists all traversal options applicable for the Linked Cells container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allLCApplicableTraversals() {
  static const std::set<TraversalOption> v {
    TraversalOption::c01, TraversalOption::c08, TraversalOption::c18, TraversalOption::sliced
#if defined(AUTOPAS_CUDA)
        ,
        TraversalOption::c01Cuda
#endif
  };
  return v;
}

/**
 * Lists all traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVCLApplicableTraversals() {
  // traversal not used but prevents usage of newton3
  static const std::set<TraversalOption> v{TraversalOption::c01};
  return v;
}

/**
 * Lists all traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVLApplicableTraversals() {
  static const std::set<TraversalOption> v{TraversalOption::verletTraversal};
  return v;
}

/**
 * Lists all traversal options applicable for the Verlet Lists Cells container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVLCApplicableTraversals() {
  static const std::set<TraversalOption> v{TraversalOption::slicedVerlet, TraversalOption::c18Verlet,
                                           TraversalOption::c01Verlet};
  return v;
}

/**
 * Lists all traversal options applicable for the given container.
 * @param container ContainerOption
 * @return set of all applicable traversal options.
 */
static inline const std::set<TraversalOption> &allApplicableTraversals(ContainerOption container) {
  switch (container) {
    case ContainerOption::linkedCells: {
      return allLCApplicableTraversals();
    }
    case ContainerOption::directSum: {
      return allDSApplicableTraversals();
    }
    case ContainerOption::verletClusterLists: {
      return allVCLApplicableTraversals();
    }
    case ContainerOption::verletLists: {
      return allVLApplicableTraversals();
    }
    case ContainerOption::verletListsCells: {
      return allVLCApplicableTraversals();
    }
  }
}

}  // namespace applicableTraversals
}  // namespace autopas