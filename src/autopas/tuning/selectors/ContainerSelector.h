/**
 * @file ContainerSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/octree/Octree.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
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
   * Constructor for the ContainerSelecor class.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff Cutoff radius to be used in this container.
   */
  ContainerSelector(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double cutoff)
      : _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff) {}

  /**
   * Container factory method.
   * @param containerChoice container to generate
   * @param containerInfo additional parameter for the container
   * @return Smartpointer to the new container
   */
  std::unique_ptr<ParticleContainerInterface<Particle_T>> generateContainer(ContainerOption containerChoice,
                                                                            ContainerSelectorInfo containerInfo);

 private:
  std::array<double, 3> _boxMin, _boxMax;
  const double _cutoff;
  ContainerSelectorInfo _currentInfo{};
};

template <class Particle_T>
std::unique_ptr<ParticleContainerInterface<Particle_T>> ContainerSelector<Particle_T>::generateContainer(
    ContainerOption containerChoice, ContainerSelectorInfo containerInfo) {
  std::unique_ptr<ParticleContainerInterface<Particle_T>> container;
  switch (containerChoice) {
    case ContainerOption::directSum: {
      container = std::make_unique<DirectSum<Particle_T>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                                          containerInfo.verletRebuildFrequency);
      break;
    }

    case ContainerOption::linkedCells: {
      container = std::make_unique<LinkedCells<Particle_T>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                                            containerInfo.verletRebuildFrequency,
                                                            containerInfo.cellSizeFactor, containerInfo.loadEstimator);
      break;
    }
    case ContainerOption::linkedCellsReferences: {
      container = std::make_unique<LinkedCellsReferences<Particle_T>>(
          _boxMin, _boxMax, _cutoff, containerInfo.verletSkin, containerInfo.verletRebuildFrequency,
          containerInfo.cellSizeFactor);
      break;
    }
    case ContainerOption::verletLists: {
      container = std::make_unique<VerletLists<Particle_T>>(
          _boxMin, _boxMax, _cutoff, containerInfo.verletSkin, containerInfo.verletRebuildFrequency,
          VerletLists<Particle_T>::BuildVerletListType::VerletSoA, containerInfo.cellSizeFactor);
      break;
    }
    case ContainerOption::verletListsCells: {
      container = std::make_unique<VerletListsCells<Particle_T, VLCAllCellsNeighborList<Particle_T>>>(
          _boxMin, _boxMax, _cutoff, containerInfo.verletSkin, containerInfo.verletRebuildFrequency,
          containerInfo.cellSizeFactor, containerInfo.loadEstimator, VerletListsCellsHelpers::VLCBuildType::soaBuild);
      break;
    }
    case ContainerOption::verletClusterLists: {
      container = std::make_unique<VerletClusterLists<Particle_T>>(
          _boxMin, _boxMax, _cutoff, containerInfo.verletSkin, containerInfo.verletRebuildFrequency,
          containerInfo.verletClusterSize, containerInfo.loadEstimator);
      break;
    }
    case ContainerOption::varVerletListsAsBuild: {
      container = std::make_unique<VarVerletLists<Particle_T, VerletNeighborListAsBuild<Particle_T>>>(
          _boxMin, _boxMax, _cutoff, containerInfo.verletSkin, containerInfo.verletRebuildFrequency,
          containerInfo.cellSizeFactor);
      break;
    }

    case ContainerOption::pairwiseVerletLists: {
      container = std::make_unique<VerletListsCells<Particle_T, VLCCellPairNeighborList<Particle_T>>>(
          _boxMin, _boxMax, _cutoff, containerInfo.verletSkin, containerInfo.verletRebuildFrequency,
          containerInfo.cellSizeFactor, containerInfo.loadEstimator, VerletListsCellsHelpers::VLCBuildType::soaBuild);
      break;
    }
    case ContainerOption::octree: {
      container =
          std::make_unique<Octree<Particle_T>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                               containerInfo.verletRebuildFrequency, containerInfo.cellSizeFactor);
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
