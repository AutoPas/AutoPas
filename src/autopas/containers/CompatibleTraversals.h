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

namespace autopas::compatibleTraversals {

/**
 * Lists all traversal options applicable for the Direct Sum container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allDSCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::ds_sequential};
  return s;
}

/**
 * Lists all traversal options applicable for the Linked Cells container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allLCCompatibleTraversals() {
  static const std::set<TraversalOption> s {
    TraversalOption::lc_c01, TraversalOption::lc_c08, TraversalOption::lc_c18, TraversalOption::lc_sliced,
        TraversalOption::lc_sliced_balanced, TraversalOption::lc_c01_combined_SoA, TraversalOption::lc_c04,
        TraversalOption::lc_c04_combined_SoA, TraversalOption::lc_c04_HCP
#if defined(AUTOPAS_CUDA)
        ,
        TraversalOption::lc_c01_cuda
#endif
  };
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Cluster Lists container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVCLCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vcl_cluster_iteration, TraversalOption::vcl_c06,
                                           TraversalOption::vcl_c01_balanced};
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVLCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vl_list_iteration};
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists Cells container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVLCCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vlc_sliced, TraversalOption::vlc_c18,
                                           TraversalOption::vlc_c01, TraversalOption::vlc_sliced_balanced};
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Cluster Cells container.
 * @return Set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVCCCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vcc_cluster_iteration_cuda};
  return s;
}

/**
 * Lists all traversal options applicable for the Var Verlet Lists As Build container.
 * @return set of all applicable traversal options.
 */
static const std::set<TraversalOption> &allVarVLAsBuildCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vvl_as_built};
  return s;
}

/**
 * Lists all traversal options applicable for the given container.
 * @param containerOption ContainerOption
 * @return set of all applicable traversal options.
 */
static inline const std::set<TraversalOption> &allCompatibleTraversals(ContainerOption containerOption) {
  switch (containerOption) {
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
    case ContainerOption::verletClusterCells: {
      return allVCCCompatibleTraversals();
    }
    case ContainerOption::varVerletListsAsBuild: {
      return allVarVLAsBuildCompatibleTraversals();
    }
  }

  autopas::utils::ExceptionHandler::exception("CompatibleTraversals: Unknown container option {}!",
                                              containerOption.to_string());

  static const std::set<TraversalOption> s{};
  return s;
}

/**
 * Lists all container options which given traversal can be applied to.
 * @param traversalOption TraversalOption
 * @return set of all compatible container options.
 */
static inline std::set<ContainerOption> allCompatibleContainers(TraversalOption traversalOption) {
  std::set<ContainerOption> result{};

  for (const auto &container : ContainerOption::getAllOptions()) {
    auto allCompatible = compatibleTraversals::allCompatibleTraversals(container);
    if (allCompatible.find(traversalOption) != allCompatible.end()) {
      result.insert(container);
    }
  }

  return result;
}

}  // namespace autopas::compatibleTraversals
