/**
 * @file StaticSelectors.h
 * @author seckler
 * @date 21.06.18
 */

#pragma once

#include <memory>

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/ReferenceLinkedCells.h"
#include "autopas/containers/verletClusterCells/VerletClusterCells.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"

namespace autopas {
/**
 * Will execute the passed function body with the static container type of container.
 * Currently, either LinkedCells, VerletLists, VerletListsCells, VerletClusterLists or DirectSum.
 *
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam FunctionType
 * @param container The container to be used.
 * @param function The function body to be executed. Has to take exactly one argument being a pointer to the container.
 * E.g: [&](auto *container){container->doSth();}  // The * is optional here. The auto is necessary!
 */
template <typename Particle, typename FunctionType>
void withStaticContainerType(std::shared_ptr<CellBasedParticleContainer<Particle>> &container, FunctionType &&function) {
  auto container_ptr = container.get();
  switch (container->getContainerType()) {
    case ContainerOption::directSum:
      function(dynamic_cast<autopas::DirectSum<Particle> *>(container_ptr));
      return;
    case ContainerOption::linkedCells:
      function(dynamic_cast<autopas::LinkedCells<Particle> *>(container_ptr));
      return;
    case ContainerOption::linkedCellsReferences:
      function(dynamic_cast<autopas::ReferenceLinkedCells<Particle> *>(container_ptr));
      return;
    case ContainerOption::verletLists:
      function(dynamic_cast<autopas::VerletLists<Particle> *>(container_ptr));
      return;
    case ContainerOption::verletListsCells:
      function(dynamic_cast<autopas::VerletListsCells<Particle> *>(container_ptr));
      return;
    case ContainerOption::verletClusterLists:
      function(dynamic_cast<autopas::VerletClusterLists<Particle> *>(container_ptr));
      return;
    case ContainerOption::verletClusterCells:
      function(dynamic_cast<autopas::VerletClusterCells<Particle> *>(container_ptr));
      return;
  }
  autopas::utils::ExceptionHandler::exception("Unknown type of container in StaticSelectorMacros.h. Type: {}",
                                              container->getContainerType());
}

}  // namespace autopas