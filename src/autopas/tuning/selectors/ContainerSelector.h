/**
 * @file ContainerSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/containers/CellBasedParticleContainer.h"
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
#include "autopas/containers/P3M/P3M_container.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/tuning/selectors/ContainerSelectorInfo.h"
#include "autopas/utils/NumParticlesEstimator.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * Selector for a particle container.
 *
 * The class is given a list of allowed container and traversal options to choose from.
 * This class selects the optimal container and delegates the choice of the optimal traversal down to this container.
 *
 * @tparam Particle_T
 * @tparam ParticleCell
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
      : _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff), _currentContainer(nullptr), _currentInfo() {}

  /**
   * Sets the container to the given option.
   * @param containerOption container to generate
   * @param containerInfo additional parameter for the container
   */
  void selectContainer(ContainerOption containerOption, ContainerSelectorInfo containerInfo);

  /**
   * Set new  boundaries, convert the container and reevaluate particle ownership.
   * @param boxMin
   * @param boxMax
   */
  void resizeBox(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
    _boxMin = boxMin;
    _boxMax = boxMax;

    _currentContainer = std::move(generateContainer(_currentContainer->getContainerType(), _currentInfo));
  }

  /**
   * Getter for the optimal container. If no container is chosen yet the first allowed is selected.
   * @return Reference to the optimal container.
   */
  inline autopas::ParticleContainerInterface<Particle_T> &getCurrentContainer();

  /**
   * Getter for the optimal container. If no container is chosen yet the first allowed is selected.
   * @return Reference to the optimal container.
   */
  inline const autopas::ParticleContainerInterface<Particle_T> &getCurrentContainer() const;

 private:
  /**
   * Container factory that also copies all particles to the new container
   * @param containerChoice container to generate
   * @param containerInfo additional parameter for the container
   * @return Smartpointer to new container
   */
  std::unique_ptr<autopas::ParticleContainerInterface<Particle_T>> generateContainer(
      ContainerOption containerChoice, ContainerSelectorInfo containerInfo);

  std::array<double, 3> _boxMin, _boxMax;
  const double _cutoff;
  std::unique_ptr<autopas::ParticleContainerInterface<Particle_T>> _currentContainer;
  ContainerSelectorInfo _currentInfo;
};

template <class Particle_T>
std::unique_ptr<autopas::ParticleContainerInterface<Particle_T>> ContainerSelector<Particle_T>::generateContainer(
    ContainerOption containerChoice, ContainerSelectorInfo containerInfo) {
  std::unique_ptr<autopas::ParticleContainerInterface<Particle_T>> container;
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
    case ContainerOption::p3m: {
      // TODO need to update default values for N_dims and cao
      //std::array<int,3> N{16,8,4};
      container = 
          std::make_unique<P3M_container<Particle_T>>(_boxMin, _boxMax, containerInfo.grid_dims, _cutoff, containerInfo.verletSkin,
                                                      containerInfo.verletRebuildFrequency,
                                                      containerInfo.cellSizeFactor, containerInfo.loadEstimator, containerInfo.cao, containerInfo.alpha);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("ContainerSelector: Container type {} is not a known type!",
                                         containerChoice.to_string());
    }
  }

  // copy particles so they do not get lost when container is switched
  if (_currentContainer != nullptr) {
    // with these assumptions slightly more space is reserved as numParticlesTotal already includes halos
    const auto numParticlesTotal = _currentContainer->size();
    const auto numParticlesHalo = autopas::utils::NumParticlesEstimator::estimateNumHalosUniform(
        numParticlesTotal, _currentContainer->getBoxMin(), _currentContainer->getBoxMax(),
        _currentContainer->getInteractionLength());

    container->reserve(numParticlesTotal, numParticlesHalo);
    for (auto particleIter = _currentContainer->begin(IteratorBehavior::ownedOrHalo); particleIter.isValid();
         ++particleIter) {
      // add particle as inner if it is owned
      if (particleIter->isOwned()) {
        container->addParticle(*particleIter);
      } else {
        container->addHaloParticle(*particleIter);
      }
    }
  }

  return container;
}

template <class Particle_T>
autopas::ParticleContainerInterface<Particle_T> &ContainerSelector<Particle_T>::getCurrentContainer() {
  if (_currentContainer == nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "ContainerSelector: getCurrentContainer() called before any container was selected!");
  }
  return *_currentContainer;
}

template <class Particle_T>
const autopas::ParticleContainerInterface<Particle_T> &ContainerSelector<Particle_T>::getCurrentContainer() const {
  if (_currentContainer == nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "ContainerSelector: getCurrentContainer() called before any container was selected!");
  }
  return *_currentContainer;
}

template <class Particle_T>
void ContainerSelector<Particle_T>::selectContainer(ContainerOption containerOption,
                                                    ContainerSelectorInfo containerInfo) {
  // Only do something if we have no container, a new type is required, or the info changed
  if (_currentContainer == nullptr or _currentContainer->getContainerType() != containerOption or
      _currentInfo != containerInfo) {
    _currentContainer = std::move(generateContainer(containerOption, containerInfo));
    _currentInfo = containerInfo;
  }
}
}  // namespace autopas
