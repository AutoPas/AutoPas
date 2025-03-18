/**
 * @file StaticContainerSelector.h
 * @author seckler
 * @date 21.06.18
 */

#pragma once

#include <memory>

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/octree/Octree.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"

namespace autopas {
/**
 * Will execute the passed function body with the static container type of container.
 *
 * @tparam ParticleT
 * @tparam ParticleCell
 * @tparam FunctionType
 * @param container The container to be used.
 * @param function The function body to be executed. Has to take exactly one argument: A reference to the container.
 * E.g: [&](auto &container){container->doSth();}  // The & is optional here. The auto is necessary!
 * @return Returns whatever function returns.
 */
template <typename ParticleT, typename FunctionType>
decltype(auto) withStaticContainerType(autopas::ParticleContainerInterface<ParticleT> &container,
                                       FunctionType &&function) {
  switch (container.getContainerType()) {
    case ContainerOption::directSum:
      return function(dynamic_cast<autopas::DirectSum<ParticleT> &>(container));
    case ContainerOption::linkedCells:
      return function(dynamic_cast<autopas::LinkedCells<ParticleT> &>(container));
    case ContainerOption::linkedCellsReferences:
      return function(dynamic_cast<autopas::LinkedCellsReferences<ParticleT> &>(container));
    case ContainerOption::verletLists:
      return function(dynamic_cast<autopas::VerletLists<ParticleT> &>(container));
    case ContainerOption::verletListsCells:
      return function(
          dynamic_cast<autopas::VerletListsCells<ParticleT, VLCAllCellsNeighborList<ParticleT>> &>(container));
    case ContainerOption::verletClusterLists:
      return function(dynamic_cast<autopas::VerletClusterLists<ParticleT> &>(container));
    case ContainerOption::pairwiseVerletLists:
      return function(
          dynamic_cast<autopas::VerletListsCells<ParticleT, VLCCellPairNeighborList<ParticleT>> &>(container));
    case ContainerOption::varVerletListsAsBuild:
      return function(
          dynamic_cast<autopas::VarVerletLists<ParticleT, VerletNeighborListAsBuild<ParticleT>> &>(container));
    case ContainerOption::octree:
      return function(dynamic_cast<autopas::Octree<ParticleT> &>(container));
  }
  autopas::utils::ExceptionHandler::exception("Unknown type of container in StaticContainerSelector.h. Type: {}",
                                              container.getContainerType());
}

}  // namespace autopas