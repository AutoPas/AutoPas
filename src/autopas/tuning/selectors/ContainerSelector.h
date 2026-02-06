/**
 * @file ContainerSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/octree/Octree.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListsSoA.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCAllCellsNeighborList.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/tuning/selectors/ContainerSelectorInfo.h"

namespace autopas {

/**
 * Selector for a particle container.
 *
 * @tparam Particle_T
 */
template <class Particle_T>
class ContainerSelector {
 public:
  /**
   * Container factory method.
   * @param containerChoice container to generate
   * @param containerInfo additional parameter for the container
   * @return Smartpointer to the new container
   */
  static std::unique_ptr<ParticleContainerInterface<Particle_T>> generateContainer(
      ContainerOption containerChoice, const ContainerSelectorInfo &containerInfo);
};

template <class Particle_T>
std::unique_ptr<ParticleContainerInterface<Particle_T>> ContainerSelector<Particle_T>::generateContainer(
    ContainerOption containerChoice, const ContainerSelectorInfo &containerInfo) {
  const auto &boxMin = containerInfo.boxMin;
  const auto &boxMax = containerInfo.boxMax;
  const auto &cutoff = containerInfo.cutoff;
  const auto &verletSkin = containerInfo.verletSkin;
  const auto &verletClusterSize = containerInfo.verletClusterSize;
  const auto &cellSizeFactor = containerInfo.cellSizeFactor;
  const auto &loadEstimator = containerInfo.loadEstimator;
  const auto &sortingThreshold = containerInfo.sortingThreshold;
  const auto &orderCellsByMortonIndex = containerInfo.orderCellsByMortonIndex;
  const auto &useOptimizedLJFunctor = containerInfo.useOptimizedLJFunctor;
  const auto &useCompactSoA = containerInfo.useCompactSoA;
  const auto &reserveVLSizes = containerInfo.reserveVLSizes;
  const auto &bucketSortParticles = containerInfo.bucketSortParticles;
  const auto &sortVerletLists = containerInfo.sortVerletLists;
  const auto &sortingFrequency = containerInfo.sortingFrequency;

  std::unique_ptr<ParticleContainerInterface<Particle_T>> container;
  switch (containerChoice) {
    case ContainerOption::directSum: {
      container = std::make_unique<DirectSum<Particle_T>>(boxMin, boxMax, cutoff, verletSkin, sortingThreshold);
      break;
    }

    case ContainerOption::linkedCells: {
      container = std::make_unique<LinkedCells<Particle_T>>(boxMin, boxMax, cutoff, verletSkin, cellSizeFactor,
                                                            sortingThreshold, loadEstimator);
      break;
    }
    case ContainerOption::linkedCellsReferences: {
      container = std::make_unique<LinkedCellsReferences<Particle_T>>(boxMin, boxMax, cutoff, verletSkin,
                                                                      cellSizeFactor, sortingThreshold);
      break;
    }
    case ContainerOption::verletLists: {
      container = std::make_unique<VerletLists<Particle_T>>(
          boxMin, boxMax, cutoff, verletSkin, VerletLists<Particle_T>::BuildVerletListType::VerletSoA, cellSizeFactor);
      container->setUseOptimizedLJFunctor(useOptimizedLJFunctor);
      break;
    }
    case ContainerOption::verletListsCells: {
      container = std::make_unique<VerletListsCells<Particle_T, VLCAllCellsNeighborList<Particle_T>>>(
          boxMin, boxMax, cutoff, verletSkin, cellSizeFactor, loadEstimator,
          VerletListsCellsHelpers::VLCBuildType::soaBuild);
      break;
    }
    case ContainerOption::verletClusterLists: {
      container = std::make_unique<VerletClusterLists<Particle_T>>(boxMin, boxMax, cutoff, verletSkin,
                                                                   verletClusterSize, loadEstimator);
      break;
    }
    case ContainerOption::varVerletListsAsBuild: {
      container = std::make_unique<VarVerletLists<Particle_T, VerletNeighborListAsBuild<Particle_T>>>(
          boxMin, boxMax, cutoff, verletSkin, cellSizeFactor);
      break;
    }

    case ContainerOption::pairwiseVerletLists: {
      container = std::make_unique<VerletListsCells<Particle_T, VLCCellPairNeighborList<Particle_T>>>(
          boxMin, boxMax, cutoff, verletSkin, cellSizeFactor, loadEstimator,
          VerletListsCellsHelpers::VLCBuildType::soaBuild);
      break;
    }
    case ContainerOption::verletListsSoA: {
      container = std::make_unique<VerletListsSoA<Particle_T>>(
      boxMin, boxMax, cutoff, verletSkin, VerletListsSoA<Particle_T>::BuildVerletListType::VerletSoA, cellSizeFactor);
      container->setUseOptimizedLJFunctor(useOptimizedLJFunctor);
      container->setUseIndexInSoAId(useCompactSoA);
      container->setReserveVLSizes(reserveVLSizes);
      container->setOrderCellsByMortonIndex(orderCellsByMortonIndex);
      container->setBucketSortParticles(bucketSortParticles);
      container->setSortVerletLists(sortVerletLists);
      container->setsortingFrequency(sortingFrequency);
      break;
    }
    case ContainerOption::octree: {
      container =
          std::make_unique<Octree<Particle_T>>(boxMin, boxMax, cutoff, verletSkin, cellSizeFactor, sortingThreshold);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("ContainerSelector: Container type {} is not a known type!",
                                         containerChoice.to_string());
    }
  }

  return container;
}
}  // namespace autopas
