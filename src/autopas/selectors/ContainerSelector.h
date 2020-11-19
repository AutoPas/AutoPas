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
#include "autopas/containers/verletClusterCells/VerletClusterCells.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/VarVerletLists.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/selectors/ContainerSelectorInfo.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * Selector for a particle container.
 *
 * The class is given a list of allowed container and traversal options to choose from.
 * This class selects the optimal container and delegates the choice of the optimal traversal down to this container.
 *
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle>
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
   * Getter for the optimal container. If no container is chosen yet the first allowed is selected.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<autopas::ParticleContainerInterface<Particle>> getCurrentContainer();

  /**
   * Getter for the optimal container. If no container is chosen yet the first allowed is selected.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<const autopas::ParticleContainerInterface<Particle>> getCurrentContainer() const;

 private:
  /**
   * Container factory that also copies all particles to the new container
   * @param containerChoice container to generate
   * @param containerInfo additional parameter for the container
   * @return smartpointer to new container
   */
  std::unique_ptr<autopas::ParticleContainerInterface<Particle>> generateContainer(ContainerOption containerChoice,
                                                                                   ContainerSelectorInfo containerInfo);

  const std::array<double, 3> _boxMin, _boxMax;
  const double _cutoff;
  std::shared_ptr<autopas::ParticleContainerInterface<Particle>> _currentContainer;
  ContainerSelectorInfo _currentInfo;
};

template <class Particle>
std::unique_ptr<autopas::ParticleContainerInterface<Particle>> ContainerSelector<Particle>::generateContainer(
    ContainerOption containerChoice, ContainerSelectorInfo containerInfo) {
  std::unique_ptr<autopas::ParticleContainerInterface<Particle>> container;
  switch (containerChoice) {
    case ContainerOption::directSum: {
      container = std::make_unique<DirectSum<Particle>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin);
      break;
    }

    case ContainerOption::linkedCells: {
      container = std::make_unique<LinkedCells<Particle>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                                          containerInfo.cellSizeFactor, containerInfo.loadEstimator);
      break;
    }
    case ContainerOption::linkedCellsReferences: {
      container = std::make_unique<LinkedCellsReferences<Particle>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                                                    containerInfo.cellSizeFactor);
      break;
    }
    case ContainerOption::verletLists: {
      container = std::make_unique<VerletLists<Particle>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                                          VerletLists<Particle>::BuildVerletListType::VerletSoA,
                                                          containerInfo.cellSizeFactor);
      break;
    }
    case ContainerOption::verletListsCells: {
      container = std::make_unique<VerletListsCells<Particle>>(_boxMin, _boxMax, _cutoff, TraversalOption::lc_c08,
                                                               containerInfo.verletSkin, containerInfo.cellSizeFactor,
                                                               containerInfo.loadEstimator);
      break;
    }
    case ContainerOption::verletClusterLists: {
      container =
          std::make_unique<VerletClusterLists<Particle>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                                         containerInfo.verletClusterSize, containerInfo.loadEstimator);
      break;
    }
    case ContainerOption::verletClusterCells: {
      container = std::make_unique<VerletClusterCells<Particle>>(_boxMin, _boxMax, _cutoff, containerInfo.verletSkin,
                                                                 containerInfo.verletClusterSize);
      break;
    }
    case ContainerOption::varVerletListsAsBuild: {
      container = std::make_unique<VarVerletLists<Particle, VerletNeighborListAsBuild<Particle>>>(
          _boxMin, _boxMax, _cutoff, containerInfo.verletSkin);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("ContainerSelector: Container type {} is not a known type!",
                                         containerChoice.to_string());
    }
  }

  // copy particles so they do not get lost when container is switched
  if (_currentContainer != nullptr) {
    for (auto particleIter = _currentContainer->begin(IteratorBehavior::haloAndOwned); particleIter.isValid();
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

template <class Particle>
std::shared_ptr<autopas::ParticleContainerInterface<Particle>> ContainerSelector<Particle>::getCurrentContainer() {
  if (_currentContainer == nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "ContainerSelector: getCurrentContainer() called before any container was selected!");
  }
  return _currentContainer;
}

template <class Particle>
std::shared_ptr<const autopas::ParticleContainerInterface<Particle>> ContainerSelector<Particle>::getCurrentContainer()
    const {
  if (_currentContainer == nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "ContainerSelector: getCurrentContainer() called before any container was selected!");
  }
  return _currentContainer;
}

template <class Particle>
void ContainerSelector<Particle>::selectContainer(ContainerOption containerOption,
                                                  ContainerSelectorInfo containerInfo) {
  // if we already have this container do nothing.
  if (_currentContainer == nullptr or _currentContainer->getContainerType() != containerOption or
      _currentInfo != containerInfo) {
    _currentContainer = std::move(generateContainer(containerOption, containerInfo));
    _currentInfo = containerInfo;
  }
}
}  // namespace autopas
