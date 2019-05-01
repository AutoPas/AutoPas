/**
 * @file ContainerSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <vector>
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/options/ContainerOption.h"
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
template <class Particle, class ParticleCell>
class ContainerSelector {
 public:
  /**
   * Constructor for the ContainerSelecor class.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff Cutoff radius to be used in this container.
   * @param cellSizeFactor Cell size factor to be used in this container.
   * @param verletSkin Length added to the cutoff for the verlet lists' skin.
   * @param verletRebuildFrequency Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   */
  ContainerSelector(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double cutoff,
                    double cellSizeFactor, double verletSkin, unsigned int verletRebuildFrequency)
      : _boxMin(boxMin),
        _boxMax(boxMax),
        _cutoff(cutoff),
        _cellSizeFactor(cellSizeFactor),
        _verletSkin(verletSkin),
        _verletRebuildFrequency(verletRebuildFrequency),
        _currentContainer(nullptr) {}

  /**
   * Sets the container to the given option.
   * @param containerOption
   */
  void selectContainer(ContainerOption containerOption);

  /**
   * Getter for the optimal container. If no container is chosen yet the first allowed is selected.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getCurrentContainer();

 private:
  /**
   * Container factory that also copies all particles to the new container
   * @param containerChoice container to generate
   * @return smartpointer to new container
   */
  std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>> generateContainer(
      ContainerOption containerChoice);

  std::array<double, 3> _boxMin, _boxMax;
  double _cutoff;
  double _cellSizeFactor;
  double _verletSkin;
  unsigned int _verletRebuildFrequency;
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> _currentContainer;
};

template <class Particle, class ParticleCell>
std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::generateContainer(ContainerOption containerChoice) {
  std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>> container;

  switch (containerChoice) {
    case directSum: {
      container = std::make_unique<DirectSum<Particle, ParticleCell>>(_boxMin, _boxMax, _cutoff);
      break;
    }
    case linkedCells: {
      container = std::make_unique<LinkedCells<Particle, ParticleCell>>(_boxMin, _boxMax, _cutoff, _cellSizeFactor);
      break;
    }
    case verletLists: {
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container =
          std::make_unique<VerletLists<Particle>>(_boxMin, _boxMax, _cutoff, _verletSkin, _verletRebuildFrequency);
      break;
    }
    case verletListsCells: {
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container = std::make_unique<VerletListsCells<Particle>>(_boxMin, _boxMax, _cutoff, TraversalOption::c08,
                                                               _verletSkin, _verletRebuildFrequency);
      break;
    }
    case verletClusterLists: {
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container = std::make_unique<VerletClusterLists<Particle>>(_boxMin, _boxMax, _cutoff, _verletSkin,
                                                                 _verletRebuildFrequency);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("Container type {} is not a known type!",
                                         utils::StringUtils::to_string(containerChoice));
    }
  }

  // copy particles so they do not get lost when container is switched
  if (_currentContainer != nullptr) {
    for (auto particleIter = _currentContainer->begin(IteratorBehavior::haloAndOwned); particleIter.isValid();
         ++particleIter) {
      // try to add every particle as inner. If it fails try as a halo.
      try {
        container->addParticle(*particleIter);
      } catch (const autopas::utils::ExceptionHandler::AutoPasException &e) {
        container->addHaloParticle(*particleIter);
      }
    }
  }

  return container;
}

template <class Particle, class ParticleCell>
std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::getCurrentContainer() {
  if (_currentContainer == nullptr) {
    autopas::utils::ExceptionHandler::exception(
        "ContainerSelector: getCurrentContainer() called before any container was selected!");
  }
  return _currentContainer;
}

template <class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::selectContainer(ContainerOption containerOption) {
  // if we already have this container do nothing.
  if (_currentContainer == nullptr) {
    _currentContainer = generateContainer(containerOption);
  }

  if (_currentContainer and _currentContainer->getContainerType() != containerOption) {
    _currentContainer = std::move(generateContainer(containerOption));
  }
}
}  // namespace autopas
