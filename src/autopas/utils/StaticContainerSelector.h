/**
 * @file StaticContainerSelector.h
 * @author seckler
 * @date 21.06.18
 */

#pragma once

#include "autopas/containers/ParticleContainerInterface.h"
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
 * @tparam Particle
 * @tparam FunctionType
 * @param container The container to be used.
 * @param function The function body to be executed. Has to take exactly one argument: A reference to the container.
 * E.g: [&](auto &container){container->doSth();}  // The & is optional here. The auto is necessary!
 * @return Returns whatever function returns.
 */
template <typename Particle, typename FunctionType>
decltype(auto) withStaticContainerType(autopas::ParticleContainerInterface<Particle> &container,
                                       FunctionType &&function) {
  switch (container.getContainerType()) {
    case ContainerOption::directSum:
      return function(dynamic_cast<autopas::DirectSum<Particle> &>(container));
    case ContainerOption::linkedCells:
      return function(dynamic_cast<autopas::LinkedCells<Particle> &>(container));
    case ContainerOption::linkedCellsReferences:
      return function(dynamic_cast<autopas::LinkedCellsReferences<Particle> &>(container));
    case ContainerOption::verletLists:
      return function(dynamic_cast<autopas::VerletLists<Particle> &>(container));
    case ContainerOption::verletListsCells:
      return function(
          dynamic_cast<autopas::VerletListsCells<Particle, VLCAllCellsNeighborList<Particle>> &>(container));
    case ContainerOption::verletClusterLists:
      return function(dynamic_cast<autopas::VerletClusterLists<Particle> &>(container));
    case ContainerOption::pairwiseVerletLists:
      return function(
          dynamic_cast<autopas::VerletListsCells<Particle, VLCCellPairNeighborList<Particle>> &>(container));
    case ContainerOption::varVerletListsAsBuild:
      return function(
          dynamic_cast<autopas::VarVerletLists<Particle, VerletNeighborListAsBuild<Particle>> &>(container));
    case ContainerOption::octree:
      return function(dynamic_cast<autopas::Octree<Particle> &>(container));
  }
  autopas::utils::ExceptionHandler::exception("Unknown type of container in StaticContainerSelector.h. Type: {}",
                                              container.getContainerType());
}

}  // namespace autopas