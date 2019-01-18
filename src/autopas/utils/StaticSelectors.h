/**
 * @file StaticSelectorMacros.h
 * @author seckler
 * @date 21.06.18
 */

#pragma once

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"

namespace autopas {
/**
 * Will execute the passed function body with the static container type of container.
 * Currently, either LinkedCells, VerletLists, VerletListsCells, VerletClusterLists or DirectSum.
 *
 * @tparam ContainerT
 * @tparam FunctionType
 * @param containerI The container to be used.
 * @param function The function body to be executed. Has to take exactly one argument being a pointer to the container.
 * E.g: [&](auto *container){container->doSth();}  // the * is optional here.
 *
 */
template <typename ContainerT, typename FunctionType>
void withStaticContainerType(ContainerT &containerI, FunctionType &&function) {
  auto container_ptr = containerI.get();
  switch (containerI->getContainerType()) {
    case ContainerOptions::directSum:
      function(
          dynamic_cast<autopas::DirectSum<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType,
                                          typename std::remove_pointer_t<decltype(container_ptr)>::ParticleCellType> *>(
              container_ptr));
      return;
    case ContainerOptions ::linkedCells:
      function(dynamic_cast<
               autopas::LinkedCells<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType,
                                    typename std::remove_pointer_t<decltype(container_ptr)>::ParticleCellType> *>(
          container_ptr));
      return;
    case ContainerOptions ::verletLists:
      function(
          dynamic_cast<autopas::VerletLists<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType> *>(
              container_ptr));
      return;
    case ContainerOptions ::verletListsCells:
      function(dynamic_cast<
               autopas::VerletListsCells<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType> *>(
          container_ptr));
      return;
    case ContainerOptions ::verletClusterLists:
      function(dynamic_cast<
               autopas::VerletClusterLists<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType> *>(
          container_ptr));
      return;
  }
  autopas::utils::ExceptionHandler::exception("wrong type of container in StaticSelectorMacros.h");
}

}  // namespace autopas