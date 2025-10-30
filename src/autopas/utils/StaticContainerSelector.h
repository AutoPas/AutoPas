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
#include "autopas/containers/verletListsCellBased/varVerletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCAllCellsNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"

namespace autopas {
/**
 * Will execute the passed function body with the static container type of container.
 *
 * @tparam Particle_T
 * @tparam FunctionType
 * @param container The container to be used.
 * @param function The function body to be executed. Has to take exactly one argument: A reference to the container.
 * E.g: [&](auto &container){container->doSth();}  // The & is optional here. The auto is necessary!
 * @return Returns whatever function returns.
 */
template <typename Particle_T, typename FunctionType>
decltype(auto) withStaticContainerType(ParticleContainerInterface<Particle_T> &container, FunctionType &&function) {
  switch (container.getContainerType()) {
    case ContainerOption::directSum:
      return function(dynamic_cast<DirectSum<Particle_T> &>(container));
    case ContainerOption::linkedCells:
      return function(dynamic_cast<LinkedCells<Particle_T> &>(container));
    case ContainerOption::linkedCellsReferences:
      return function(dynamic_cast<LinkedCellsReferences<Particle_T> &>(container));
    case ContainerOption::verletLists:
      return function(dynamic_cast<VerletLists<Particle_T> &>(container));
    case ContainerOption::verletListsCells:
      return function(dynamic_cast<VerletListsCells<Particle_T, VLCAllCellsNeighborList<Particle_T>> &>(container));
    case ContainerOption::verletClusterLists:
      return function(dynamic_cast<VerletClusterLists<Particle_T> &>(container));
    case ContainerOption::pairwiseVerletLists:
      return function(dynamic_cast<VerletListsCells<Particle_T, VLCCellPairNeighborList<Particle_T>> &>(container));
    case ContainerOption::varVerletListsAsBuild:
      return function(dynamic_cast<VarVerletLists<Particle_T, VerletNeighborListAsBuild<Particle_T>> &>(container));
    case ContainerOption::octree:
      return function(dynamic_cast<Octree<Particle_T> &>(container));
  }
  utils::ExceptionHandler::exception("Unknown type of container in StaticContainerSelector.h. Type: {}",
                                     container.getContainerType());
}

}  // namespace autopas