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
#include "autopas/containers/verletListsCellBased/verletListsCells/PseudoVerletLists.h"

namespace autopas {
/**
 * Will execute the passed function body with the static container type of container.
 *
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam FunctionType
 * @param container The container to be used.
 * @param function The function body to be executed. Has to take exactly one argument being a pointer to the container.
 * E.g: [&](auto *container){container->doSth();}  // The * is optional here. The auto is necessary!
 * @return Returns whatever function returns.
 */
template <typename Particle, typename FunctionType>
decltype(auto) withStaticContainerType(const std::shared_ptr<ParticleContainerInterface<Particle>> &container,
                                       FunctionType &&function) {
  auto containerPtr = container.get();
  switch (container->getContainerType()) {
    case ContainerOption::directSum:
      return function(dynamic_cast<autopas::DirectSum<Particle> *>(containerPtr));
    case ContainerOption::linkedCells:
      return function(dynamic_cast<autopas::LinkedCells<Particle> *>(containerPtr));
    case ContainerOption::linkedCellsReferences:
      return function(dynamic_cast<autopas::LinkedCellsReferences<Particle> *>(containerPtr));
    case ContainerOption::verletLists:
      return function(dynamic_cast<autopas::VerletLists<Particle> *>(containerPtr));
    case ContainerOption::verletListsCells:
      return function(
          dynamic_cast<autopas::VerletListsCells<Particle, VLCAllCellsNeighborList<Particle>> *>(containerPtr));
    case ContainerOption::verletClusterLists:
      return function(dynamic_cast<autopas::VerletClusterLists<Particle> *>(containerPtr));
    case ContainerOption::pairwiseVerletLists:
      return function(
          dynamic_cast<autopas::VerletListsCells<Particle, VLCCellPairNeighborList<Particle>> *>(containerPtr));
    case ContainerOption::varVerletListsAsBuild:
      return function(
          dynamic_cast<autopas::VarVerletLists<Particle, VerletNeighborListAsBuild<Particle>> *>(containerPtr));
    case ContainerOption::octree:
      return function(dynamic_cast<autopas::Octree<Particle> *>(containerPtr));
    case ContainerOption::pseudoVerletLists:
      return function(dynamic_cast<autopas::PseudoVerletLists<Particle, PVLCellPairNeighborList<Particle>> *>(containerPtr));
  }
  autopas::utils::ExceptionHandler::exception("Unknown type of container in StaticContainerSelector.h. Type: {}",
                                              container->getContainerType());
}

}  // namespace autopas