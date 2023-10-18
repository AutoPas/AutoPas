/**
 * @file CompatibleTraversals.h
 * @author F. Gratl
 * @date 5/31/19
 */

#pragma once

#include <set>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"

namespace autopas::compatibleTraversals {

/**
 * Lists all traversal options applicable for the Direct Sum container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allDSCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::ds_sequential};
  return s;
}

/**
 * Lists all 3-body traversal options applicable for the Direct Sum container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allDSCompatibleTraversals3B() {
  static const std::set<TraversalOption> s{TraversalOption::ds_sequential_3b};
  return s;
}

/**
 * Lists all traversal options applicable for the Linked Cells container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allLCCompatibleTraversals() {
  static const std::set<TraversalOption> s{
      TraversalOption::lc_c01,
      TraversalOption::lc_c08,
      TraversalOption::lc_c18,
      TraversalOption::lc_sliced,
      TraversalOption::lc_sliced_balanced,
      TraversalOption::lc_c01_combined_SoA,
      TraversalOption::lc_c04,
      TraversalOption::lc_c04_combined_SoA,
      TraversalOption::lc_c04_HCP,
      TraversalOption::lc_sliced_c02,
  };
  return s;
}

/**
 * Lists all 3-body traversal options applicable for the Linked Cells container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allLCCompatibleTraversals3B() {
  static const std::set<TraversalOption> s{TraversalOption::lc_c01_3b};
  return s;
}

/**
 **
 * Lists all traversal options applicable for the Reference Linked Cells container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allRLCCompatibleTraversals() {
  return allLCCompatibleTraversals();
}

/**
 * Lists all traversal options applicable for the Verlet Cluster Lists container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVCLCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vcl_cluster_iteration, TraversalOption::vcl_c06,
                                           TraversalOption::vcl_c01_balanced,      TraversalOption::vcl_sliced,
                                           TraversalOption::vcl_sliced_balanced,   TraversalOption::vcl_sliced_c02};
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vl_list_iteration};
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists Cells container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLCCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vlc_sliced, TraversalOption::vlc_c18,
                                           TraversalOption::vlc_c01, TraversalOption::vlc_sliced_c02,
                                           TraversalOption::vlc_sliced_balanced};
  return s;
}

/**
 * Lists all traversal options applicable for the Var Verlet Lists As Build container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVarVLAsBuildCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vvl_as_built};
  return s;
}

/**
 * Lists all traversal options applicable for the Pairwise Verlet Lists container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLPCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::vlp_sliced,
                                           TraversalOption::vlp_c18,
                                           TraversalOption::vlp_c01,
                                           TraversalOption::vlp_sliced_c02,
                                           TraversalOption::vlp_sliced_balanced,
                                           TraversalOption::vlp_c08};
  return s;
}

/**
 * Lists all traversal options applicable for the Octree container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allOTCompatibleTraversals() {
  static const std::set<TraversalOption> s{TraversalOption::ot_c18, TraversalOption::ot_c01};
  return s;
}

/**
 * Provides a set of all traversals that only support Newton3 mode disabled.
 * @return
 */
[[maybe_unused]] static std::set<TraversalOption> allTraversalsSupportingOnlyNewton3Disabled() {
  return {TraversalOption::lc_c01,
          TraversalOption::lc_c01_combined_SoA,
          TraversalOption::ot_c01,
          TraversalOption::vcl_c01_balanced,
          TraversalOption::vcl_cluster_iteration,
          TraversalOption::vl_list_iteration,
          TraversalOption::vlc_c01,
          TraversalOption::vlp_c01,
          TraversalOption::lc_c01_3b};
};
/**
 * Provides a set of all traversals that only support Newton3 mode enabled.
 * @return
 */
[[maybe_unused]] static std::set<TraversalOption> allTraversalsSupportingOnlyNewton3Enabled() {
  return {
      TraversalOption::ot_c18,
  };
};
/**
 * Provides a set of all traversals that only support DataLayout AoS.
 * @return
 */
[[maybe_unused]] static std::set<TraversalOption> allTraversalsSupportingOnlyAoS() { return {}; };
/**
 * Provides a set of all traversals that only support DataLayout SoA.
 * @return
 */
[[maybe_unused]] static std::set<TraversalOption> allTraversalsSupportingOnlySoA() {
  return {
      TraversalOption::lc_c01_combined_SoA,
      TraversalOption::lc_c04_combined_SoA,
  };
};

/**
 * Lists all traversal options applicable for the given container.
 * @param containerOption ContainerOption
 * @param interactionTypeOption Interaction type for which compatible traversals are collected
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allCompatibleTraversals(
    ContainerOption containerOption, const InteractionTypeOption interactionTypeOption) {
  switch (interactionTypeOption) {
    // Check compatible pairwise traversals
    case InteractionTypeOption::pairwise: {
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
        case ContainerOption::varVerletListsAsBuild: {
          return allVarVLAsBuildCompatibleTraversals();
        }
        case ContainerOption::linkedCellsReferences: {
          return allRLCCompatibleTraversals();
        }
        case ContainerOption::pairwiseVerletLists: {
          return allVLPCompatibleTraversals();
        }
        case ContainerOption::octree: {
          return allOTCompatibleTraversals();
        }
      }
    }
    // Check compatible 3-body traversals
    case InteractionTypeOption::threeBody: {
      switch (containerOption) {
        case ContainerOption::directSum: {
          return allDSCompatibleTraversals3B();
        }
        case ContainerOption::linkedCells: {
          return allLCCompatibleTraversals3B();
        }
        default: {
          static const std::set<TraversalOption> s{};
          return s;
        }
      }
    }
  }

  autopas::utils::ExceptionHandler::exception(
      "CompatibleTraversals: Unknown container option {} or unsupported interaction type {}!",
      containerOption.to_string(), interactionTypeOption.to_string());

  static const std::set<TraversalOption> s{};
  return s;
}

/**
 * Lists all container options which given traversal can be applied to.
 * @param traversalOption TraversalOption
 * @return set of all compatible container options.
 */
[[maybe_unused]] static std::set<ContainerOption> allCompatibleContainers(TraversalOption traversalOption) {
  std::set<ContainerOption> result{};

  for (const auto &container : ContainerOption::getAllOptions()) {
    for (const auto &interactionType : InteractionTypeOption::getAllOptions()) {
      auto allCompatible = compatibleTraversals::allCompatibleTraversals(container, interactionType);
      if (allCompatible.find(traversalOption) != allCompatible.end()) {
        result.insert(container);
      }
    }
  }

  return result;
}

}  // namespace autopas::compatibleTraversals
