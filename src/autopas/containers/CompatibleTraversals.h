/**
 * @file CompatibleTraversals.h
 * @author F. Gratl
 * @date 5/31/19
 */

#pragma once

#include <set>
#include "autopas/options/ContainerOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {
namespace compatibleTraversals {

/**
 * Lists all traversal options applicable for the Direct Sum container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allDSCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::directSumTraversal};
  return s;
}

/**
 * Lists all traversal options applicable for the Linked Cells container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allLCCompatibleTraversals() {
  static const std::set<TraversalOption> s {
    TraversalOption::c01, TraversalOption::c08, TraversalOption::c18, TraversalOption::sliced,
        TraversalOption::c01CombinedSoA, TraversalOption::c04, TraversalOption::c04SoA
#if defined(AUTOPAS_CUDA)
        ,
        TraversalOption::c01Cuda
#endif
  };
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVCLCompatibleTraversals() {
  // traversal not used but prevents usage of newton3
  static const std::set<TraversalOption> s{TraversalOption::c01};
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVLCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::verletTraversal};
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists Cells container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVLCCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::slicedVerlet, TraversalOption::c18Verlet,
                                           TraversalOption::c01Verlet};
  return s;
}

/**
 * Lists all traversal options applicable for the given container.
 * @param container ContainerOption
 * @return set of all applicable traversal options.
 */
static inline const std::set<TraversalOption> &allCompatibleTraversals(ContainerOption container) {
  switch (container) {
    case ContainerOption::linkedCells: {
      return allLCCompatibleTraversals();
    }
    case ContainerOption::directSum: {
      return allDSCompatibleTraversals();
    }
    case ContainerOption::verletClusterLists: {
      return allVCLCompatibleTraversals();
    }
    case ContainerOption::verletLists: {
      return allVLCompatibleTraversals();
    }
    case ContainerOption::verletListsCells: {
      return allVLCCompatibleTraversals();
    }
  }

  autopas::utils::ExceptionHandler::exception("OptionSelector: Unknown selector strategy {}!",
                                              autopas::utils::StringUtils::to_string(container));

  static const std::set<TraversalOption> s{};
  return s;
}

}  // namespace compatibleTraversals
}  // namespace autopas