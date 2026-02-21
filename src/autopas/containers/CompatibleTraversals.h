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
 * Anonymous namespace for internal helper functions.
 */
namespace {
/**
 * Helper function to filter all traversal options for a given prefix and interaction type.
 * @param prefix
 * @param interactionType
 * @return Set of all options that match the prefix and interaction type.
 */
std::set<TraversalOption> filterAllOptions(const std::string &prefix, const InteractionTypeOption &interactionType) {
  const auto allOpts = TraversalOption::getAllOptionsOf(interactionType);
  std::set<TraversalOption> retSet;
  // If the lambda condition holds (=prefix matches) copy the option in the return set.
  std::copy_if(allOpts.begin(), allOpts.end(), std::inserter(retSet, retSet.begin()),
               [&](const auto &option) { return option.to_string().substr(0, prefix.length()) == prefix; });
  return retSet;
}
}  // namespace

/**
 * Lists all traversal options applicable for the Direct Sum container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allDSCompatibleTraversals() {
  static const auto s = filterAllOptions("ds_", InteractionTypeOption::pairwise);
  return s;
}

/**
 * Lists all triwise traversal options applicable for the Direct Sum container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allDSCompatibleTraversals3B() {
  static const auto s = filterAllOptions("ds_", InteractionTypeOption::triwise);
  return s;
}

/**
 * Lists all traversal options applicable for the Linked Cells container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allLCCompatibleTraversals() {
  static const auto s = filterAllOptions("lc_", InteractionTypeOption::pairwise);
  return s;
}

/**
 * Lists all triwise traversal options applicable for the Linked Cells container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allLCCompatibleTraversals3B() {
  static const auto s = filterAllOptions("lc_", InteractionTypeOption::triwise);
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
  static const auto s = filterAllOptions("vcl_", InteractionTypeOption::pairwise);
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLCompatibleTraversals() {
  static const auto s = filterAllOptions("vl_", InteractionTypeOption::pairwise);
  return s;
}

/**
 * Lists all triwise traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLCompatibleTraversals3B() {
  static const auto s = filterAllOptions("vl_", InteractionTypeOption::triwise);
  return s;
}

/**
 * Lists all traversal options applicable for the Verlet Lists Cells container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLCCompatibleTraversals() {
  static const auto s = filterAllOptions("vlc_", InteractionTypeOption::pairwise);
  return s;
}

/**
 * Lists all triwise traversal options applicable for the Verlet Lists container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLCCompatibleTraversals3B() {
  static const auto s = filterAllOptions("vlc_", InteractionTypeOption::triwise);
  return s;
}

/**
 * Lists all traversal options applicable for the Var Verlet Lists As Build container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVarVLAsBuildCompatibleTraversals() {
  static const auto s = filterAllOptions("vvl_", InteractionTypeOption::pairwise);
  return s;
}

/**
 * Lists all traversal options applicable for the Pairwise Verlet Lists container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allVLPCompatibleTraversals() {
  static const auto s = filterAllOptions("vlp_", InteractionTypeOption::pairwise);
  return s;
}

/**
 * Lists all traversal options applicable for the Octree container.
 * @return set of all applicable traversal options.
 */
[[maybe_unused]] static const std::set<TraversalOption> &allOTCompatibleTraversals() {
  static const auto s = filterAllOptions("ot_", InteractionTypeOption::pairwise);
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
          // TraversalOption::vl_list_iteration,
          TraversalOption::vl_list_intersection_sorted_3b,
          TraversalOption::vl_list_intersection_hashing_3b,
          TraversalOption::vl_pair_list_iteration_3b,
          TraversalOption::vlc_c01,
          TraversalOption::vlp_c01};
}
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
[[maybe_unused]] static std::set<TraversalOption> allTraversalsSupportingOnlyAoS() {
  return {
      TraversalOption::vl_list_intersection_sorted_3b,
  };
};
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
    // Check compatible triwise traversals
    case InteractionTypeOption::triwise: {
      switch (containerOption) {
        case ContainerOption::directSum: {
          return allDSCompatibleTraversals3B();
        }
        case ContainerOption::linkedCells: {
          return allLCCompatibleTraversals3B();
        }
        case ContainerOption::verletLists: {
          return allVLCompatibleTraversals3B();
        }
        case ContainerOption::verletListsCells: {
          return allVLCCompatibleTraversals3B();
        }
        default: {
          static const std::set<TraversalOption> s{};
          return s;
        }
      }
    }
    case InteractionTypeOption::all: {
      utils::ExceptionHandler::exception("CompatibleTraversals: Interaction type option {} is not supported!",
                                         interactionTypeOption.to_string());
    }
  }

  utils::ExceptionHandler::exception(
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
    for (const auto &interactionType : InteractionTypeOption::getMostOptions()) {
      auto allCompatible = compatibleTraversals::allCompatibleTraversals(container, interactionType);
      if (allCompatible.find(traversalOption) != allCompatible.end()) {
        result.insert(container);
      }
    }
  }

  return result;
}

}  // namespace autopas::compatibleTraversals
