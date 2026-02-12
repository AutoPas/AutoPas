/**
 * @file CompatibleCellSizeFactors.h
 * @date 12.02.2026
 * @author S. J. Newcome
 */

#pragma once

#include <set>

#include "autopas/options/ContainerOption.h"


namespace autopas::compatibleCSFs {
/**
 * Lists all containers only supporting CSF 1 (or ignore CSF).
 * @return set of all containers that only support CSF 1 (or ignore CSF)
 */
[[maybe_unused]] static const std::set<ContainerOption> allContainersOnlySupportingCSF1() {
  return {ContainerOption::octree, ContainerOption::verletClusterLists, ContainerOption::directSum};
}

/**
 * Lists all containers that support CSF < 1.0. Note, it is possible that a container supports this but one or more
 * traversals of that container do not. These should be filtered out in <traversal>::isApplicableToDomain().
 * @return A list of all containers that support CSF < 1.0.
 */
[[maybe_unused]] static const std::set<ContainerOption> allContainersSupportingSub1CSF() {
  return {ContainerOption::linkedCells, ContainerOption::linkedCellsReferences};
}

}  // namespace autopas::compatibleCSFs